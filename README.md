# GFoil

Viscous-inviscid aerofoil solver with boundary-layer analysis, trailing-edge noise prediction (Amiet/WPS), and CoDiPack-based automatic differentiation.

---

## Build

Requires: CMake ≥ 3.15, C++17 compiler, Eigen 3, pybind11, Python 3.8.

```bash
# Configure (once)
cmake -B build . \
  -DPYBIND11_PYTHON_VERSION=3.8 \
  -DPYTHON_EXECUTABLE=$(pyenv which python3.8)

# Build all targets
cmake --build build -j$(nproc)
```

This produces:
| Target | Location | Purpose |
|--------|----------|---------|
| `GFoil_fwd_codi` | `build/` | Standalone forward solver binary |
| `GFoil_AD` | `build/` | Standalone AD (gradient) solver binary |
| `gfoil_cpp.so` | `GFoil/` | pybind11 module for the Python API |

Build a single target:
```bash
cmake --build build --target gfoil_cpp   # pybind11 module only
```

---

## Python API

```python
import numpy as np
from GFoil import fwd_run, grad_run
from GFoil.inputs import Aerofoil, Acoustics, OperatingConds

# Load airfoil coordinates (any .dat file with x, y columns)
coords = np.loadtxt("n0012_sharp.dat")
foil = Aerofoil(xcoords=coords[:, 0], ycoords=coords[:, 1])

op = OperatingConds(alpha=2.0, Re=2e6, nCrit=9.0)
ac = Acoustics(observerXYZ=[0.0, 3.0, 0.5])

result = fwd_run(foil, op, ac)
print(result.CL, result.CD, result.OASPL)

# Gradients (dCL/dy, dCD/dy, dOASPL/dy, plus alpha scalars)
grads = grad_run(result, foil, op, ac)
```

---

## Regression test

```bash
# After a fresh build (golden files already committed)
python3 tests/regression_test.py --test

# Rebuild then test
python3 tests/regression_test.py --build --test
```

Golden files are in `tests/golden/`. The canonical test input (`tests/input.json`) is copied to the repo root automatically before each run.

To regenerate golden files after an intentional physics change:
```bash
python3 tests/regression_test.py --create-golden
git add tests/golden/ tests/input.json
git commit -m "update golden reference after <describe change>"
```

---

## Airfoil data files

The root directory contains `.dat` files in two-column `x y` format (one coordinate pair per line, 300–400 points, starting at the trailing edge). These are passed directly to `Aerofoil(xcoords=..., ycoords=...)` after loading with `numpy.loadtxt`.

---

## Environment variables

| Variable | Effect |
|----------|--------|
| `GFOIL_DEBUG=1` | Enable per-iteration diagnostics (residual, omega, ilam, SparseLU NaN reports) |
