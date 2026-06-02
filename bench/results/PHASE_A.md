# Part A — low-Re acoustic NaN fixed at the source

## Result

| metric | pre-A | post-A |
|---|---|---|
| full cold-converged (valid CL/CD/CM **and** OASPL) | 1184/1375 = 86.1 % | **1203/1375 = 87.5 %** |
| `acoustic_nan` failures | 19 | **0** |
| regressions on previously-converging cases | — | **0** |

All 19 `acoustic_nan` cases now return `converged=true` with a finite OASPL. The
other failure buckets are unchanged (no_convergence 50, transition_front_oscillation
74, nan_lock 24, diverged 24). Golden regression **10/10 bit-identical** (forward
and AD). Noise code (Amiet/WPS formulas) untouched — only input conditioning.

## Diagnosis (Step A0, committed separately)

`GFOIL_DEBUG` traces of all 19 cases showed two distinct mechanisms, both at
Re=2e5 (mostly-laminar TE with at most a thin near-separated turbulent layer):

1. **13 cases — `a` overflows to `inf` in `calc_WPS_Rozenburg`.** `tauWall`
   degenerates to ~4e-5 (Cf ~ 1e-5, near separation), driving the Clauser
   parameter `beta_c = (theta/tauWall)·dpdx` to **90–920**, far outside the
   Rozenberg calibration. The amplitude exponent `A1 = 3.7 + 1.5·beta_c` then
   overflows `pow(base, A1)`. `Delta`, `Rt`, `SS` all stayed finite.
2. **6 cases — `integral = 0` → `10·log10(0) = -inf`.** Both surfaces are fully
   laminar at the TE (`tauMax = 0`, so `calc_WPS` is never called), i.e. there is
   genuinely **no TBL-TE noise source**.

## Fix (Step A1, this commit) — `src/include/sound.hpp`

**Physical input floors in `calc_WPS`** (applied for every model — roz/goo/lee/kam/tno):
- `Ue ≥ 1e-6` (existing pattern).
- `theta ≥ 1e-12`; `deltaStar ≥ 1.05·theta` (a turbulent BL has H > 1; removes the
  `Delta = delta/deltaStar → ∞` and `SS = Ue/(tauMax²·deltaStar) → ∞` paths);
  `delta ≥ deltaStar`.
- `tauWall ≥ max(Cf_min·q, theta·|dpdx|/beta_max)` with `q = ½ρUe²`, `Cf_min = 1e-4`
  (well below any attached-turbulent Cf), `beta_max = 50` (the edge of the
  empirical APG calibration). The second term is the key one: it floors the wall
  shear at exactly the value that keeps the Clauser parameter within the model's
  validity range, i.e. "the thinnest resolvable attached turbulent layer."
- `tauMax ≥ tauWall`.

These are smooth `std::max` floors and, verified under `GFOIL_DEBUG`, **never bind
for a normal attached TE BL** — at Re=1e6 the same NACA 0018 case runs with
`beta_c = 1.7 / 10.8` (uncapped) and unfloored `tauWall`, so all Re≥1e6 and the
Re=2e6 golden are bit-identical. The floor only engages in the degenerate low-Re
near-separation regime, where the result is a bounded estimate (the empirical WPS
model is out of calibration there — treat low-Re acoustic numbers as
floor-limited / low-confidence).

**Acoustic silence floor** for the zero-source case: branch on the passive value
so a source-carrying integral takes the *exact original* `10·log10(integral/pref²)`
(bit-identical tape + adjoint — an `std::max` here perturbs `dOASPL/dy` at ~1e-7),
while a zero/negative integral returns a constant `10·log10(1e-30) ≈ -300 dB`. A
fully-laminar TE genuinely radiates no modelled TBL-TE noise, so −300 dB is the
honest "no source" limit; its OASPL adjoint is correctly **zero** (no source ⇒ no
acoustic shape sensitivity).

## Validation

- OASPL is **monotonic in Re** (no non-physical reversal): e.g. NACA 0018 α=−3°
  nCrit=12 reads −2.2 dB (Re=2e5, floor-limited) → 48.3 (1e6) → 86.0 (5e6);
  AG10 α=0° reads −300 (laminar TE, no source) → 46.0 (1e6, turbulent TE) → 80.1
  (5e6). The laminar→turbulent step is a real regime change, not an artefact.
- **OASPL adjoint differentiates cleanly** for floored cases (finite `dOASPL/dy`,
  finite `dOASPL/dα`) and is exactly zero for the silence-floor cases.
- `run_forward.cpp`'s `acoustic_nan` label is retained as a backstop; it is now
  rarely (never, on this grid) hit.

Artifacts: `sweep_postA.jsonl`, `sweep_postA.log`; diagnostic `bench/diag_acoustic_nan.py`.
