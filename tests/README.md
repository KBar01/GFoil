# GFoil Regression Tests

## Overview

The suite drives the `GFoil.gfoil_cpp` pybind11 module directly (no standalone
binaries, no `restart.json`) and runs five groups of checks:

1. **Free-transition golden case** — `tests/input.json` (NACA 0012 sharp-TE
   coordinates, alpha=2°, Re=2e6, nCrit=5, model `roz`, observer (0,3,0.5),
   span 2): forward scalars (CL, CD, CM, OASPL) plus AD alpha scalars and all
   three gradient arrays. Tolerance 1e-8 relative.
2. **Forced-transition golden case** — same conditions with transition forced
   at x/c = 0.1 on both surfaces. Exercises the taped-`xift` adjoint path
   (CHANGELOG "Forced-transition adjoint fix"), which previously had no
   regression coverage. Same outputs and tolerance; golden files carry the
   `forced_` prefix.
3. **Coarse-input golden case** — analytic open-TE NACA 0012 with only 101
   nodes (`naca0012_analytic(n_half=50)`), free-transition conditions of
   group 1. Exercises the runtime input-geometry length (`Nin` removed, see
   CHANGELOG "Variable-length input geometry"); gradient arrays are length
   101. Golden files carry the `coarse101_` prefix.
4. **AD-vs-FD spot check** — on the coarse case: central finite differences
   (h=1e-5, forward `rtol=1e-11`) at 5 y-nodes × {CL, CD, OASPL} against the
   AD gradients, tolerance 2e-3 (typical agreement ~1e-5). Nodes deliberately
   avoid the transition-sensitive x≈0.2 region where the response is kinked
   on the ±1e-5 scale (pre-existing physics; see the CHANGELOG entry). Plus a
   `noise_run` smoke check driven by the coarse case's TE BL states.
5. **Windowed-Amiet anchors** — two forward-only OASPL reference values from
   CHANGELOG "Windowed Amiet TE sample / Reference values": scalar
   `TEsample=0.98` → 63.18598 dB and window `[0.95,0.99]` → 63.36702 dB,
   tolerance 1e-6 relative. **The anchor config is fully embedded in
   `regression_test.py`** and uses an *analytic open-TE NACA 0012* (cosine
   spacing, 301 points), model `kam`, observer (1,0,1), span 3, rtol 1e-6 —
   NOT the `input.json` coordinates (those give 63.05499/63.24251; the foil
   generation is part of the recorded config).

48 checks total. Exit code is **0** if all pass, **1** otherwise.

---

## Other suites

These are separate scripts, same conventions (plain script, PASS/FAIL, exit 0
iff all pass). They are **not** run by `regression_test.py`.

| Script | Covers |
|---|---|
| `tests/amiet_kernel_test.py --test` | The general oblique-gust Amiet kernel (29 checks): algebra identities, cut behaviour and the near-cutoff bridge, golden `golden/amiet_kernel_scalars.json` |
| `tests/rotor_noise_test.py --test` | The rotor TE-noise wrapper (44 checks), golden `golden/rotor_noise_scalars.json` |
| `tests/midspan_reduction_check.py` | One-off before/after gate (below) |

### The mid-span reduction gate

`midspan_reduction_check.py` is not a golden test — it compares two *builds*.
The general `_vec` kernel must reduce to the mid-span kernel it replaced, at
x2 = 0, to rtol ≤ 1e-13. Run it whenever the `_vec` kernel changes:

```bash
git stash push src/include/newAmiet.hpp            # or check out the old kernel
cmake --build build --target gfoil_cpp -j8
python3 tests/midspan_reduction_check.py --dump /tmp/ref.json
git stash pop
cmake --build build --target gfoil_cpp -j8
python3 tests/midspan_reduction_check.py --check /tmp/ref.json
```

Last measured: worst relative error 8.6e-15 over 42 spectra × 64 frequencies,
45% of values bit-identical. See CHANGELOG "General oblique-gust Amiet kernel"
for why the remainder is not bit-identical (it is `Radiation_integral2_general`'s
`D`, `Y²` and `coeffI`, algebraically identical at κ̄ = μ̄ but differently
associated).

---

## Quick reference

```bash
# Regenerate goldens — ONLY after an intentional numerics change (see below)
python3 tests/regression_test.py --create-golden

# After every refactoring step
python3 tests/regression_test.py --test

# Rebuild then test
python3 tests/regression_test.py --build --test

# Or via the shell wrapper (from the repo root)
bash tests/run_regression.sh --test
```

---

## Tolerance

The default relative tolerance for golden comparisons is **1e-8**. A value
passes if

```
rel_err = |v_current - v_golden| / |v_golden|   (or absolute if |v_golden| < 1e-14)
rel_err < tol
```

Override with `--tol 1e-6` if needed. The anchor checks use their own fixed
1e-6 tolerance (the recorded values have 7 significant figures).

## Updating the golden files

If a **legitimate algorithmic change** (not a pure refactor) intentionally
alters the numerics, regenerate on the updated build, and record the
justification in CHANGELOG.md:

```bash
python3 tests/regression_test.py --create-golden
git add tests/golden/
git commit -m "update golden reference after <describe change>"
```

Do **not** update golden files to paper over a regression. Forward scalars are
expected to be *bit-identical* (rel_err exactly 0) across pure refactors and
build-flag changes; gradient arrays may shift at the ~1e-7 level for changes
that reorder reverse-mode accumulation (see CHANGELOG for precedents).

---

## Files

| File | Purpose |
|------|---------|
| `tests/regression_test.py` | Main test script (golden cases + embedded anchor configs) |
| `tests/run_regression.sh` | Shell wrapper for convenience |
| `tests/input.json` | Canonical golden-case input (NACA 0012 sharp-TE coords) |
| `tests/golden/fwd_scalars.json` | Golden forward scalars (free transition) |
| `tests/golden/ad_scalars.json` | Golden AD alpha scalars (free) |
| `tests/golden/ad_gradients.json` | Golden AD gradient arrays (free) |
| `tests/golden/forced_fwd_scalars.json` | Golden forward scalars (forced transition) |
| `tests/golden/forced_ad_scalars.json` | Golden AD alpha scalars (forced) |
| `tests/golden/forced_ad_gradients.json` | Golden AD gradient arrays (forced) |
| `tests/golden/coarse101_fwd_scalars.json` | Golden forward scalars (101-node input) |
| `tests/golden/coarse101_ad_scalars.json` | Golden AD alpha scalars (101-node input) |
| `tests/golden/coarse101_ad_gradients.json` | Golden AD gradient arrays (101-node, length 101) |
| `tests/amiet_kernel_test.py` | General oblique-gust Amiet kernel suite |
| `tests/golden/amiet_kernel_scalars.json` | Golden \|I\| sweeps (both branches + bridge, rtol 1e-10) |
| `tests/rotor_noise_test.py` | Rotor TE-noise wrapper suite |
| `tests/golden/rotor_noise_scalars.json` | Golden rotor OASPL / pinned Spp (rtol 1e-10) |
| `tests/midspan_reduction_check.py` | Mid-span reduction gate (before/after builds) |
