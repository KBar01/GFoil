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

- **`Aerofoil(xcoords, ycoords, chord=1.0, span=...)`** — any number of
  coordinate points `n >= 10` is accepted (no upper limit); the geometry is
  spline-fit and re-panelled internally onto a fixed 200-node distribution, so
  input coarseness affects only how well the spline captures the true shape.
  The trailing edge must be at `x = 1.0` on both surfaces. Gradient arrays
  returned by `grad_run` have length `n`.
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

## Rotor noise post-processing (`rotor_noise.py`)

Predicts the time-averaged far-field **broadband trailing-edge noise** of a
rotor by wrapping the aerofoil solver's acoustics in the corrected
Schlinker–Amiet rotating-blade procedure of Sinayoko, Kingan & Agarwal (2013),
*Proc. R. Soc. A* 469:20130065. It is pure Python post-processing — no C++, no
aero solve of its own, never on the AD/optimisation path.

You supply the blade as radial **strips**, each with its own trailing-edge
boundary-layer state from your own `fwd_run` at that strip's sectional
conditions. The wrapper does the rotor kinematics and acoustics.

### The procedure

Per (strip, azimuth, observer), in the hub-fixed frame (rotor turns about `+z`
at `Omega`; axial inflow `Uz` runs in `-z`):

1. **Emission time** — solve `c0*Te = |xo - xe - M_FO*c0*Te|` for the retarded
   time. (Generalises the paper's far-field form: the source position stays
   finite, so moderate observer distances are valid.)
2. **Source positions** — convected `xc` and *present* `xp` (Eq. 4.1).
3. **Doppler ratio** — from the convected source-to-observer direction
   (Eq. 4.10); the source frequencies are `omega' = omega/doppler`.
4. **Blade frame** — rotate the observer about `xp` into the chord-aligned
   section frame (Eq. 4.9), and evaluate the fixed-aerofoil PSD there at
   `omega'`.
5. **Assemble** — weight by `(omega'/omega)**2`, average over azimuth, sum
   strips incoherently, multiply by the blade count (Eq. 4.12).

The wall-pressure spectrum is a function of *source* frequency, so each call is
handed the Doppler-shifted grid directly.

### Example

```python
import numpy as np
from GFoil import (fwd_run, Aerofoil, OperatingConds, Acoustics,
                   RotorConfig, RotorStrip, rotor_noise_run,
                   strip_relative_speed)

cfg = RotorConfig(Omega=120.0, Uz=15.0, B=3)     # rad/s, m/s, blades

strips = []
for radius, chord, pitch_deg, chi_deg in [(1.2, 0.16, 12.0, 4.0),
                                          (1.6, 0.14, 9.0, 3.5),
                                          (2.0, 0.12, 6.0, 3.0)]:
    # One forward solve per strip, at ITS sectional conditions.
    U_rel = strip_relative_speed(cfg, radius)
    foil = Aerofoil(x, y, chord=chord, span=0.4)          # x, y as in Quick start
    op = OperatingConds(alpha=chi_deg, Re=U_rel * chord / cfg.nu, nCrit=9.0)
    r = fwd_run(foil, op, Acoustics(observerXYZ=np.array([0.0, 0.0, 1.0])),
                verbose=True)
    if not r.converged:
        raise SystemExit(f"strip r={radius} did not converge: {r.failure_mode}")

    strips.append(RotorStrip(
        radius=radius, dr=0.4, chord=chord, pitch_rad=np.deg2rad(pitch_deg),
        BL_top=r.verbose_data.BL_top,      # pass the TE BL states straight through
        BL_bot=r.verbose_data.BL_bot,
        wps_model="kam",
    ))

freqs = np.logspace(np.log10(500), np.log10(15000), 60)   # observer-frame [Hz]
obs = np.array([[0.0, 40.0, 30.0], [0.0, 0.0, 50.0]])     # hub frame [m]

res = rotor_noise_run(cfg, strips, freqs, obs)
print(res.OASPL_perObs)                       # (nObs,)  dB re 20 uPa
print(res.Spp.shape, res.Spp_dB.shape)        # (nObs, N) linear Pa^2/(rad/s), dB/Hz
print(res.diagnostics["spanwise_neglect_ratio"])   # obliquity indicator
```

A strip can instead carry `custom_WPS_func=f`, where
`f(omega_src) -> (WPS_upper, WPS_lower)` in linear `Pa^2/(rad/s)`, evaluated at
the Doppler-shifted source frequencies. Pair it with `Ue_custom`: the Amiet
stage still needs an edge velocity (`U_c = 0.7*Ue`) even with a custom spectrum.

### Validity regime

`kC > 1`, `omega >> Omega`, and spanwise correlation length `l_S << radius`. All
three are reported in `RotorNoiseResult.diagnostics` and warn when marginal;
none ever fails a run.

### Status: off-axis directivity works (phase 2, July 2026)

The phase-1 mid-span-kernel limitation is **resolved**. `noise_run`'s Amiet
stage now uses the general three-dimensional oblique-gust formulation of Roger &
Moreau (2005): it forms `S0 = sqrt(x1^2 + beta^2*(x2^2 + x3^2))`, selects the
Eq. 18 gust `K2_bar = k_bar*x2/S0`, branches supercritical/subcritical on
`xi = beta*|x2|/S0` with a near-cutoff regularisation, and applies the
spanwise-wavenumber-corrected Corcos length. Phase 1 put the rotor plane 21.3 dB
*above* the axis — an inverted dipole; the axis is now 12.5 dB above the plane,
as the edge-dipole physics requires. Gated by `tests/rotor_noise_test.py` group
10. See CHANGELOG, "General oblique-gust Amiet kernel (phase 2)".

`diagnostics["spanwise_neglect_ratio"]` is still reported but no longer measures
an error — it is now an obliquity indicator and a monotone proxy for `xi`.

### Remaining limitations — read before trusting a level

- **Eq. 18 gust selection.** The kernel takes the paper's large-aspect-ratio
  limit, where the spanwise wavenumber integral collapses to a delta selecting
  one gust per (observer, frequency), rather than a finite-span sinc over many.
  Valid for `L/(2b) >> 1`; strips of near-unit aspect ratio are outside it.
- **Axial inflow only** — no cross-flow or shaft angle.
- **Kernel `c0` is fixed at 340 m/s.** `RotorConfig.c0` drives the wrapper's
  kinematics but cannot reach the kernel; a mismatch > 1 m/s warns.
- **No A-weighting.**
- **`fwd_run` and `noise_run` diverge off mid-span.** The taped forward/AD path
  (`calc_OASPL`) is still the mid-span kernel; only `noise_run` is general. They
  agree to 8.6e-15 for a mid-span observer. Off mid-span, `noise_run` is the one
  to trust. Optimisation observers are mid-span, so the taped path was left
  byte-identical rather than re-taped for no gradient benefit.

Tests: `python3 tests/rotor_noise_test.py --test` (44 checks; the key gates are
the static anchor — with `Omega = 0` the wrapper must reduce to a bare
`noise_run` call bit-for-bit — and the axis-vs-plane directivity ordering) and
`python3 tests/amiet_kernel_test.py --test` (29 checks on the kernel itself).

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
