# GFoil

Viscous–inviscid aerofoil solver with boundary-layer analysis, trailing-edge
noise prediction (Amiet / wall-pressure-spectrum models), and CoDiPack-based
adjoint automatic differentiation. The numerical core is C++; you drive it from
Python through a pybind11 extension module.

---

## Requirements

- Python ≥ 3.8 with `numpy`
- A C++17 compiler (GCC ≥ 9 or Clang ≥ 10)
- CMake ≥ 3.14 (Ninja optional)

You do **not** need to install the C++ libraries yourself: CMake's FetchContent
downloads Eigen, CoDiPack, nlohmann/json and pybind11 at build time. **An
internet connection is required on the first build** (the dependencies are then
cached by CMake for subsequent builds).

---

## Install

```bash
git clone <repo-url> GFoil
cd GFoil
pip install .
```

The first install takes a few minutes while CMake fetches the C++ dependencies
and compiles the extension. For development (editable install that picks up
Python-side changes without recompiling):

```bash
pip install --no-build-isolation -e .
```

`--no-build-isolation` lets scikit-build-core reuse the existing CMake build
directory for editable mode.

---

## Quick start

`fwd_run` runs the forward solve (aerodynamics + trailing-edge noise); pass its
result to `grad_run` for gradients. The snippet is self-contained — it builds an
analytic NACA 0012 — and is directly runnable.

```python
import numpy as np
from GFoil import fwd_run, grad_run, Aerofoil, OperatingConds, Acoustics

# Coordinates: arrays x, y, chord normalised to 1.0, ordered trailing-edge ->
# lower surface -> leading edge -> upper surface -> trailing edge. An XFOIL .dat
# file loads directly with:  x, y = np.loadtxt("foil.dat").T
n = 150
beta = np.linspace(0.0, np.pi, n + 1)
xs = 0.5 * (1 + np.cos(beta))
yt = lambda x: 5 * 0.12 * (0.2969*np.sqrt(x) - 0.1260*x
                           - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
x = np.concatenate([xs, xs[::-1][1:]])
y = np.concatenate([-yt(xs), yt(xs[::-1][1:])])

foil      = Aerofoil(x, y, chord=1.0, span=2.0)
operating = OperatingConds(alpha=3.0, Re=2.0e6, Ma=0.0, nCrit=9.0)
acoustics = Acoustics(observerXYZ=np.array([0.0, 0.0, 3.0]), model="kam")

result = fwd_run(foil, operating, acoustics)
if not result.converged:
    raise SystemExit(f"forward solve did not converge: {result.failure_mode}")

print(f"CL={result.CL:.4f}  CD={result.CD:.5f}  "
      f"CM={result.CM:.4f}  OASPL={result.OASPL:.2f} dB")

# Adjoint gradients: alpha sensitivities + full per-coordinate sensitivity arrays.
grads = grad_run(result, foil, operating, acoustics)
print(f"dCL/dalpha={grads.dCL_dalpha:.5f}  dOASPL/dalpha={grads.dOASPL_dalpha:.4f}")
print(f"dCL/dy: length {len(grads.dCL_dy)} (one entry per coordinate)")
```

### Inputs

- **`Aerofoil(xcoords, ycoords, chord=1.0, span=...)`** — coordinates are
  re-panelled internally; the trailing edge must be at `x = 1.0` on both
  surfaces.
- **`OperatingConds(alpha=, Re=, Ma=0.0, nCrit=9.0, rtol=1e-6)`** — `alpha` in
  degrees; `rtol` is the RMS convergence tolerance (tighten to ≤ 1e-10 for
  AD-vs-finite-difference gradient checks).
- **`Acoustics(observerXYZ=, model="kam", TESampleLoc=0.98, ...)`** —
  `observerXYZ` is `(3,)` or `(N, 3)` in the freestream-aligned frame with origin
  at quarter-chord; `model` is one of `roz`/`goo`/`lee`/`kam`/`tno`;
  `TESampleLoc` is a scalar `x/c` or a `[x_lo, x_hi]` averaging window.

### Outputs

- `fwd_run` → `FwdResult`: `.converged`, `.CL`, `.CD`, `.CM`, `.OASPL`,
  `.failure_mode` (and `.verbose_data` when called with `verbose=True`).
- `grad_run` → `GradResult`: `.dCL_dalpha`, `.dCD_dalpha`, `.dOASPL_dalpha` and
  the coordinate-sensitivity arrays `.dCL_dy`, `.dCD_dy`, `.dOASPL_dy`.

---

## Developers

- Regression suite (forward scalars + AD scalars/arrays, free & forced
  transition, windowed-Amiet anchors):

  ```bash
  python3 tests/regression_test.py --test          # compare against golden files
  python3 tests/regression_test.py --build --test  # rebuild, then test
  ```

  Golden files live in `tests/golden/`; regenerate with `--create-golden` only
  after an intentional physics change.

- `bench/` holds the cold-start convergence sweep and the baseline contract
  (`bench/results/REPORT.md`) that gates solver changes. The geometry library it
  samples (`Smoothed_TEfixed_linear/`) is a large external dataset kept out of
  the repo.

- Architecture, CoDiPack rules and known limitations live in `CLAUDE.md`; the
  chronological change record is in `CHANGELOG.md`.
