# GFoil

Viscous-inviscid aerofoil solver with boundary-layer analysis, trailing-edge noise prediction (Amiet/WPS), and CoDiPack-based automatic differentiation.

---

## Requirements

- Python ≥ 3.8
- C++17 compiler (GCC ≥ 9 or Clang ≥ 10)
- CMake ≥ 3.14
- ninja (recommended, not required)
- numpy

Eigen, CoDiPack, nlohmann/json, and pybind11 are downloaded automatically by CMake during the build — you do not need to install them.

---

## Installation

```bash
pip install .
```

The first install takes a few minutes: CMake downloads the C++ dependencies and compiles the extension module.

For development (editable install, picks up Python changes without reinstalling):

```bash
pip install --no-build-isolation -e .
```

`--no-build-isolation` is required because scikit-build-core needs access to the already-built CMake cache when working in editable mode.

---

## Quick start

```python
import numpy as np
from GFoil import fwd_run, Aerofoil, Acoustics, OperatingConds

coords = np.loadtxt("your_aerofoil.dat")  # two columns: x, y; chord normalised to 1.0

foil      = Aerofoil(coords[:, 0], coords[:, 1], chord=1.0, span=2.0)
operating = OperatingConds(alpha=3.0, Re=2e6, Ma=0.0, nCrit=9.0)
acoustics = Acoustics(observerXYZ=np.array([0.0, 0.0, 3.0]), model='kam')

result = fwd_run(foil, operating, acoustics)

if result.converged:
    print(f"CL={result.CL:.4f}  CD={result.CD:.6f}  OASPL={result.OASPL:.2f} dB")
else:
    print(f"Did not converge: {result.failure_mode}")
```

---

## Gradients

`grad_run` computes dCL/dα, dCD/dα, dOASPL/dα and the full dCL/dy, dCD/dy, dOASPL/dy coordinate-sensitivity arrays via adjoint automatic differentiation (CoDiPack). Pass the `FwdResult` from a converged forward solve:

```python
from GFoil import grad_run

grads = grad_run(result, foil, operating, acoustics)
# grads.dCL_dalpha, grads.dCD_dalpha, grads.dOASPL_dalpha
# grads.dCL_dy, grads.dCD_dy, grads.dOASPL_dy  (arrays, length = number of coordinates)
```

---

## Coordinate file format

Coordinates are two columns (x, y) with chord normalised to 1.0 and the trailing edge at x = 1.0. Both the upper and lower surface trailing-edge points must be present (x = 1.0 for each). XFOIL-format `.dat` files work directly with `numpy.loadtxt`.

---

## Regression test

```bash
python3 tests/regression_test.py --test        # test against committed golden files
python3 tests/regression_test.py --build --test # rebuild then test
```

Golden files are in `tests/golden/`. To regenerate after an intentional physics change:

```bash
python3 tests/regression_test.py --create-golden
git add tests/golden/ tests/input.json
git commit -m "update golden reference after <describe change>"
```
