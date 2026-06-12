# GFoil Regression Tests

## Overview

The suite drives the `GFoil.gfoil_cpp` pybind11 module directly (no standalone
binaries, no `restart.json`) and runs three groups of checks:

1. **Free-transition golden case** — `tests/input.json` (NACA 0012 sharp-TE
   coordinates, alpha=2°, Re=2e6, nCrit=5, model `roz`, observer (0,3,0.5),
   span 2): forward scalars (CL, CD, CM, OASPL) plus AD alpha scalars and all
   three gradient arrays. Tolerance 1e-8 relative.
2. **Forced-transition golden case** — same conditions with transition forced
   at x/c = 0.1 on both surfaces. Exercises the taped-`xift` adjoint path
   (CHANGELOG "Forced-transition adjoint fix"), which previously had no
   regression coverage. Same outputs and tolerance; golden files carry the
   `forced_` prefix.
3. **Windowed-Amiet anchors** — two forward-only OASPL reference values from
   CHANGELOG "Windowed Amiet TE sample / Reference values": scalar
   `TEsample=0.98` → 63.18598 dB and window `[0.95,0.99]` → 63.36702 dB,
   tolerance 1e-6 relative. **The anchor config is fully embedded in
   `regression_test.py`** and uses an *analytic open-TE NACA 0012* (cosine
   spacing, 301 points), model `kam`, observer (1,0,1), span 3, rtol 1e-6 —
   NOT the `input.json` coordinates (those give 63.05499/63.24251; the foil
   generation is part of the recorded config).

22 checks total. Exit code is **0** if all pass, **1** otherwise.

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
