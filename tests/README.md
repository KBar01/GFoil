# GFoil Regression Tests

## Overview

The regression suite runs both the forward solver (`GFoil_fwd_codi`) and the AD
solver (`GFoil_AD`) against `input.json`, then compares their outputs to a set of
golden reference files stored in `tests/golden/`.

---

## Quick reference

```bash
# First-time setup — run once on known-good code before any refactoring
python3 tests/regression_test.py --create-golden

# After every refactoring step
python3 tests/regression_test.py --test

# Rebuild then test
python3 tests/regression_test.py --build --test

# Or via the shell wrapper (from the repo root)
bash tests/run_regression.sh --test
```

Exit code is **0** if all checks pass, **1** if any fail.

---

## Workflow

### 1. First-time setup

Before making any code changes, run `--create-golden` on a known-good build.
This executes both solvers and saves their outputs as the reference:

```
tests/golden/fwd_scalars.json   — CL, CD, CM, OASPL
tests/golden/ad_scalars.json    — d cl/d alpha, d cd/d alpha, d OASPL/d alpha
tests/golden/ad_gradients.json  — d cl/d ycoords, d cd/d ycoords, d OASPL/d ycoords
```

**Commit the golden files to the repository** so every collaborator tests against
the same reference.

### 2. Iterative refactoring

After each templating or refactoring step, run:

```bash
python3 tests/regression_test.py --test
```

If the script exits **1**, a value drifted beyond the tolerance — revert the
last change and investigate before continuing.

### 3. Tolerance

The default relative tolerance is **1e-8**.  A value `v_current` passes if:

```
rel_err = |v_current - v_golden| / |v_golden|   (or absolute if |v_golden| < 1e-14)
rel_err < tol
```

Override with `--tol 1e-6` if needed.

### 4. Updating the golden files

If a **legitimate algorithmic change** (not a pure refactor) intentionally alters
the numerics, regenerate the golden files on the updated build:

```bash
python3 tests/regression_test.py --create-golden
git add tests/golden/
git commit -m "update golden reference after <describe change>"
```

Do **not** update golden files to paper over a regression.

---

## Solver dependency

The forward solver must run before the AD solver because the AD solver reads
`restart.json`, which the forward solver writes when `"returnData": 1` is set in
`input.json`.  The test script always runs them in the correct order.

Both binaries must be run from the repo root (where `input.json` lives).
`regression_test.py` handles this automatically via `cwd=REPO_ROOT`.

## Exit codes

Both GFoil binaries exit with code **1** even on a successful run (this appears
to be a convergence-status indicator in the solver, not an error).  The test
script therefore does **not** use the exit code as a success signal.  Instead it
checks that the expected output file (`out.json` / `ad_gradients.json`) was
created or updated after the run.  A true failure leaves the output file absent
or unchanged.

---

## Files

| File | Purpose |
|------|---------|
| `tests/regression_test.py` | Main test script |
| `tests/run_regression.sh` | Shell wrapper for convenience |
| `tests/golden/fwd_scalars.json` | Golden forward-solver scalars |
| `tests/golden/ad_scalars.json` | Golden AD scalar gradients |
| `tests/golden/ad_gradients.json` | Golden AD array gradients |
