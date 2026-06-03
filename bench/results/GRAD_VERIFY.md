# AD-vs-FD gradient verification — root-cause diagnosis & fix

**Date:** 2026-06-03 · **Branch:** `grad-verify-investigation`
**Case:** NACA0012-family SVD foil (`n0012_sharp.dat` baseline, `Smoothed_TEfixed_linear/`
library, 10 SVD modes), chord 0.3, span 1.5, Re 2e6, Ma 0, nCrit 6, α 3°, model `kam`,
observer (0,0,3), TESampleLoc 0.98.
**Drivers:** `bench/h1_diag.py`, `bench/h1_stepsweep.py`, `bench/grad_sweep.py`,
`bench/plot_before_after.py`. FD = central difference, AD = two-pass adjoint
(`partialOutputspartialInputs` + `partialRpartialx`).

---

## TL;DR

| Symptom | Hypothesis | Verdict |
|---|---|---|
| 1. Forced transition AD≠FD | **H1** (xift detached from tape) | **CONFIRMED & FIXED** |
| 2. Higher modes → larger rel error | **H1** (constant missing term / smaller gradient) | **CONFIRMED & FIXED** |
| 3. "Cliff" in error-vs-h | H2 (gated branch) | **Reclassified:** it is the flat **H1 error floor**, not a discontinuity cliff. H2 ruled out for this case. |
| (α column check) | H3 (radiation geometry) | **Ruled out** — α agrees to 0.000% before *and* after. |

The single fix (make `xift` a taped `Real`, measured against the **same**
`isol_final.distFromStag` the residual stations use) drives **every** forced-transition
gradient — CL, CD, OASPL, all 10 modes and α — to the FD noise floor, with the
regression suite remaining **bit-identical** (free-transition golden untouched).

---

## H1 — forced-transition station `xift` was detached from the AD tape (CONFIRMED)

### Mechanism
The residual `R` depends on the forced-transition arc-length position `xift` (it sets
`xt`, which weights the transition-station interpolation and enters `residual_station`
through `xt − x1`, `x2 − xt`; see `residuals_shared.hpp:502–565`). `xift` is a
*continuous* function of geometry (and, through the stagnation point, of α) via
`foil.x` and `distFromStag`.

In `ADfuncs.hpp::partialRpartialx` the entire `xift` interpolation was computed with
`.getValue()` (all `double`), and `Vsol_t::xift` / `Param_t::xift` were stored as
`double`. So `d(xift)/dy ≡ 0` and `d(xift)/dα ≡ 0` on the tape, and the adjoint term
`−λᵀ ∂R/∂x` was missing `∂R/∂xift · ∂xift/∂x`. Pass 1
(`partialOutputspartialInputs`) freezes `U` and never recomputes the transition path,
so the term was absent from **both** passes — the forced-transition dependence was
silently missing entirely.

### Confirming experiment 1 — free vs forced (`bench/h1_diag.py`, h=1e-5)
Free transition agrees well; forced is badly wrong, and CD (most transition-sensitive)
is wrecked:

```
                         OASPL%      CL%        CD%
free   (worst mode)        3.7       0.04       1.1
forced mode 4              4.4      23.5      394.6
forced mode 9             21.3       2.2      992.7
```

### Confirming experiment 2 — step-size shape (`bench/h1_stepsweep.py`)
Forced-transition error is a **flat floor, independent of h** (mode 4 ≈ 4.43%, mode 7 ≈
19.5%, mode 9 ≈ 21.4% across h ∈ [1e-3, 1e-7], then round-off noise). A constant
h-independent offset is the signature of a **missing constant term**, not truncation
(which falls ∝h²) and not a discontinuity cliff (which would be high-then-step-down).

### Why symptom 2 (mode-number dependence) follows from H1
The missing term `∂R/∂xift·∂xift/∂w_i` is an *absolute* error of similar size across
modes; the *relative* error is therefore larger where the gradient magnitude is smaller,
i.e. at the higher (lower-energy) SVD modes. This is exactly what the data show.

---

## The α subtlety — why a naïve fix broke `dOASPL/dα` (and how it was resolved)

A first fix that taped `xift` against `isol_pre.distFromStag` (the **inviscid**
stagnation solution) *improved every geometry mode* but **regressed the α column**
(OASPL α 0.000% → 3.4%, CD α 0.000% → 11%).

Root cause: `distFromStag[i] = |foil.s[i] − stagArcLocation|`
(`solver_funcs.hpp:325`) moves rigidly with the stagnation point, which shifts strongly
with α. The residual depends on the transition position **relative** to its stations
(`xt − x1`), both measured from the same stag point, so an α-driven stag shift must
**cancel** between `xift` and the stations. The residual stations are built from
`isol_final.distFromStag` (post `stagpoint_move_AD`, viscous-corrected), but
`isol_pre.distFromStag` has a *different* `d(stagArc)/dα`. Mixing the two broke the
cancellation and injected a spurious `dxift/dα`.

The forward solver never has this problem: `coupled.cpp` runs
`stagpoint_move → update_transition → build_glob_RV`, so `xift` and the stations share
the **same** post-move `distFromStag` and the cancellation holds (this is why FD — and
pre-fix AD — both gave α = 0.000%). The fix mirrors the forward ordering: compute
`xift` from `isol_final.distFromStag` *after* `stagpoint_move_AD`.

---

## H2 — step-size "cliff" / gated branches (RULED OUT for this case)

Symptom 3 ("abrupt drop in error at a particular step size") is **not** a crossed
discontinuity here. The forced-transition before-fix curve is a **flat floor** at the
H1 missing-term level for all usable h (e.g. CD mode 4 = 394.58% from h=1e-3 down to
h≈1e-7), then degrades into round-off noise. After the H1 fix the same curves are clean
truncation **V**s bottoming at ≈0.

Decisive evidence that the H2 candidate branches (the −300 dB acoustic floor and
`mean_power` floor in `sound.hpp`; the `calc_WPS` `std::max` input floors; the
`Radiation_integral2` `1e-10` clips) are **not active/binding** for this configuration:
**after the H1 fix the OASPL gradient agrees with FD to 0.000% for every mode.** Had any
of those branches been binding and kinking the spectrum, an irreducible OASPL error
would survive the H1 fix. It does not. (Consistent with OASPL ≈ 75 dB ≫ the −300 dB
floor, and an attached high-Re TE BL where the `calc_WPS` floors are documented not to
bind.) Any residual *large-h* artifact that may appear for other cases is the integer
panel-bracket / turbulent-node selection in the transition logic — legitimately
non-differentiable; AD is locally correct there and FD is the contaminated reference.

---

## H3 — radiation-geometry α dependence (RULED OUT)

`alpha_rad` reaches `calc_OASPL` as the taped `Real alpha = (alphad/180)·π` in pass 1
(`ADfuncs.hpp:112`), so `cos_a/sin_a/x_loc/z_loc` carry `d/dα`. `partialRpartialx` does
not call `calc_OASPL` (the OASPL radiation-geometry dependence is an explicit partial,
owned by pass 1) — correct by construction. Empirically, the **α row agrees to 0.000%
both before and after** the H1 fix, so the radiation-direction geometry contribution is
complete. H3 is not implicated.

---

## The fix (minimal, tape-correct)

1. `src/include/data_structs_shared.hpp` — `Vsol_t::xift` and `Param_t::xift`:
   `double → Real` (taped). Forward structs in `data_structs.h` stay `double` (the
   forward solver needs only the value; bit-identical).
2. `srcAD/include/ADfuncs.hpp::partialRpartialx` — moved the `xift` computation to
   **after** `stagpoint_move_AD` and compute the interpolation in `Real` against
   `isol_final.distFromStag`. The integer **panel bracket** selection stays passive
   (`.getValue()` comparisons — legitimately non-differentiable); only the continuous
   in-bracket interpolation `xift = xi_prev + (xi_curr−xi_prev)·frac` (and `xft_abs`,
   `frac`) is taped.
3. `srcAD/include/main_func.hpp:89` — `param.xift = … : Real(0.0)` (ternary type).
   `residuals_shared.hpp:527` `static_cast<Real>(param.xift)` is now an identity and
   preserves the tape automatically.

**Guardrail note (Param_t Real field).** CLAUDE.md warns that a `Real` field with a
default initializer in `Param_t<Real>` can shift AD gradients. Empirically it does
**not** here: the regression (free transition) is bit-identical, because `Real xift =
0.0` is a passive constant and only becomes tape-active when assigned from
`isol_final.distFromStag` in forced-transition cases (which the golden does not
exercise).

---

## Results

`python3 tests/regression_test.py --build --test` → **10/10 PASS, bit-identical**
(forward CL/CD/CM/OASPL and all AD scalars/arrays, rel_err 0.000e+00) before *and*
after the fix.

Forced-transition AD-vs-FD, worst-case (min-over-h) relative error %
(`bench/plot_before_after.py`; plots `GRAD_VERIFY_{CL,CD,OASPL}.png`):

```
   row | CL before  CL after | CD before  CD after | OASPL before  OASPL after
 mode1 |     0.021     0.000 |     0.229     0.000 |       0.152       0.000
 mode2 |     3.078     0.000 |     7.586     0.000 |       3.594       0.000
 mode3 |     0.215     0.000 |   144.947     0.000 |       2.438       0.000
 mode4 |     9.067     0.000 |   391.543     0.000 |       2.713       0.000
 mode5 |     0.879     0.000 |    21.659     0.000 |       1.808       0.000
 mode6 |     0.014     0.000 |    16.119     0.000 |       0.367       0.000
 mode7 |     0.527     0.000 |    61.058     0.000 |       8.747       0.000
 mode8 |     0.250     0.000 |   169.142     0.000 |       0.233       0.000
 mode9 |     0.547     0.000 |   975.078     0.000 |       5.453       0.000
mode10 |     0.072     0.000 |     0.278     0.000 |       0.098       0.000
 alpha |     0.000     0.000 |     0.000     0.000 |       0.000       0.000
```

Every forced-transition gradient collapses to the FD noise floor; α is preserved at
0.000%. CD was hit hardest pre-fix (drag is the most transition-sensitive output).
