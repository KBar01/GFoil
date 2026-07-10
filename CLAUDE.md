### Entry point
`calc_OASPL<Real, WriteJSON>()` in `sound.hpp`:
- Builds log-spaced frequency array (200–20000 Hz, Nsound=250 points)
- Calls `calc_WPS<Real>()` dispatch for upper and lower surface WPS
- Calls `TE_noise_outer<Real>()` for far-field PSD
- Integrates to get OASPL in dB (ref 20μPa)
- WriteJSON=true only in forward build (JSON output compiled out in AD via
  `if constexpr`)

### WPS models (src/include/WPSmodels.hpp)
All template<typename Real>. Dispatch via string key in calc_WPS():
- `"roz"` — Rozenburg
- `"goo"` — Goody (2004)
- `"lee"` — Lee
- `"kam"` — Kamruzzaman
- `"tno"` — TNO

`calc_WPS` floors its physical inputs before dispatch (deltaStar≥1.05·theta,
delta≥deltaStar, tauWall≥max(Cf_min·q, theta·|dpdx|/beta_max) with Cf_min=1e-4,
beta_max=50) to keep degenerate low-Re / near-separated TE BL states inside the
empirical models' validity (otherwise beta_c→O(100s) overflows the Rozenberg
amplitude). Floors never bind for attached TE BLs (Re≥1e6, golden bit-identical).
`calc_OASPL` floors a zero acoustic source (fully-laminar TE) to a finite −300 dB
instead of log10(0)=−inf. See `bench/results/PHASE_A.md`.

### newAmiet.hpp — key functions
- `errFunc<Real>`             — erf via Faddeeva; CoDi derivatives via
                                StatementPushHelper. CoDi types ONLY.
- `Fresnel_int_conj<Real>`    — computes E*(x); E(x) = conj(E*(x)) by negating imag
- `Radiation_integral1<Real>` — R&M Eq. 13, primary TE scattering
- `Radiation_integral2<Real>` — R&M Eq. 14, back-scattering correction
- `Radiation_integral_total<Real>` — frequency loop over Nsound points
- `TE_noise_outer<Real>`      — top-level noise function

### TE_noise_outer signature (rho, nu, Ky removed in May 2026 audit)
```cpp
void TE_noise_outer(
    Real M, Real U, Real x, Real y, Real z,
    Real b,          // semi-chord = c/2
    Real c,          // full chord
    Real span,
    Real c0,         // speed of sound = 340 m/s
    const Real omega[Nsound],
    const Real Ue_bot, const Real Ue_top,
    Real (&WPS_lower)[Nsound], Real (&WPS_upper)[Nsound],
    Real (&farfieldSpectra)[Nsound])
```

### Call site in sound.hpp (line ~120)
```cpp
TE_noise_outer<Real>(c, Uinf, X, Y, Z, chordScale/2.0, chordScale,
                     S, 340.0, omega,
                     edgeVel_bot, edgeVel_top,
                     WPSLower, WPSUpper, farfieldSpectra);
// c=M, Uinf=U, X=x, Y=y, Z=z — 15 args total, exact match
```

---

## Real Type and Include Rules (CRITICAL)

| File | Defines |
|---|---|
| `src/include/real_type.h` | `Real = codi::RealReverse`, all macros incl. Nsound |
| `srcAD/include/real_type.hpp` | same macros, no Real typedef (uses template) |

**`Nin` no longer exists** (July 2026): input geometry length is a runtime
`int nIn` threaded from the bindings through `runCode`/the AD drivers to
`make_panels`/`spline_curvature`, with hard floor `NinMin = 10` (defined in
both real_type headers; mirrored as `N_MIN_INPUT_NODES` in `inputs.py`).
Gradient arrays returned by `grad_run` have length `nIn`. Internal
discretisation (`Nfine`, `Ncoords`, `RVdimension`, …) is unchanged and fixed.

**Never include both in the same TU** — `norm2` redefinition results.

**Nsound** must be defined before any noise header. Satisfied by include order
in both TUs — do NOT add `#include "real_type.h"` inside newAmiet.hpp.

### Include paths (CMakeLists.txt)
- `GFoil_fwd_codi`: `src/include/`, `src/noise_includes/`
- `GFoil_AD`: `srcAD/include/`, `src/include/`, `srcAD/noise_includes/`,
  `src/noise_includes/`

---

## CoDi Rules (CRITICAL)

| Do | Don't |
|---|---|
| `std::abs` | `std::fabs` — not specialized for CoDi |
| `std::sqrt(a*a + b*b)` | `std::hypot(a,b)` — not specialized |
| `std::sin`, `std::cos`, `std::pow`, `std::sqrt` | `using std::foo` at header scope |
| Fully-qualified `std::` names in headers | `.getValue()` outside CoDi context |

`StatementPushHelper` requires an active tape. `errFunc` uses `.getValue()` and
`.getGradient()` — valid only when Real is a CoDi AD type, never with double.

---

## CoDi Type Rules (EXTENDED)

### Never mix double and Real
Even for quantities that are not design variables (e.g. observer
coordinates, constants, loop indices), do NOT cast to double or pass
as double* in any function that is instantiated with a CoDi Real type.
Mixing double arithmetic with Real arithmetic can silently detach
computations from the CoDi tape, producing wrong gradients with no
compile error.

The only place getValue() is legitimately used is:
  - Inside errFunc() in newAmiet.hpp, where StatementPushHelper
    manually registers the erf derivative (this is the intended
    CoDi external function pattern).
  - In WriteJSON blocks (if constexpr (WriteJSON)) where we are
    extracting passive values for output only, not feeding back
    into any Real computation.

Everywhere else: keep everything as Real.

### Observer coordinates
Observer locations (X, Y, Z) are not design variables but must still
be typed as Real throughout calc_OASPL, TE_noise_outer, and all
intermediate functions. The CoDi tape simply records zero gradient
contribution from these — which is correct behaviour, not a problem.

### Passing arrays of Real
When passing observer locations or similar non-design Real arrays, use:
  const Real* obsX, const Real* obsY, const Real* obsZ
Never as double* or std::vector<double>.

### Do NOT add Real fields to Param_t<Real> in data_structs_shared.hpp
`Param_t<Real>` is instantiated inside CoDi active tape regions. Any `Real`
field with a default initializer executes a tape assignment at construction,
creating a spurious entry and silently corrupting AD gradients (confirmed:
adding `ncrithyst` to `Param_t<Real>` shifted gradient arrays by ~6×10⁻⁷).
Rule: fields needed only by the forward solver belong in the non-template
`Param` struct in `data_structs.h`, not in `Param_t`.

Exception: `Param_t::xift` / `Vsol_t::xift` ARE `Real` (taped) — deliberately, so
the forced-transition station carries `d(xift)/d(geometry)` into the adjoint (see
the forced-transition adjoint fix in CHANGELOG.md and `bench/results/GRAD_VERIFY.md`).
This is safe because the default `Real xift = 0.0` is a passive constant until it is
assigned from a taped quantity (`isol_final.distFromStag`) in forced-transition
cases only; the free-transition golden is bit-identical. Do NOT revert these to
`double` — that re-detaches the forced-transition gradient. The matching forward
structs in `data_structs.h` stay `double` (value-only, no tape).

---

## Shared solver template headers

`src/include/solver_funcs.hpp` — included by both build targets.
Contains duck-typed template implementations; each fwd *.cpp delegates
via a one-line non-template wrapper. main_func.hpp includes it directly.

`src/include/panel_funcs.hpp` — already shared; `inviscid_velocity` and
`dvelocity_dgamma` live here as `template<Real,FoilT>`.

### Dead / Stale Files
- `src/noise_includes/amiet.h` — deleted (May 2026); was not included anywhere.

---

## Testing / Regression

Golden-output regression test covers CL, CD, CM, OASPL (forward) and
all three gradient arrays + alpha scalars (AD).

```
# First-time setup (golden files committed; only needed after physics changes)
python3 tests/regression_test.py --create-golden

# Run after every refactoring step
python3 tests/regression_test.py --test

# Rebuild then test
python3 tests/regression_test.py --build --test
```

Golden files: `tests/golden/fwd_scalars.json`, `ad_scalars.json`, `ad_gradients.json`
Test input:   `tests/input.json` (NACA 0012, alpha=2°, ncrit=5)
Tolerance: 1e-8 relative.

The test drives the pybind11 module (`GFoil.gfoil_cpp`) directly — no standalone
binaries involved. The standalone `GFoil_fwd_codi` and `GFoil_AD` executables are
build artefacts only (main.cpp is a stub).

---

## Known Limitations

Full write-ups with root cause analysis and rejected approaches are in CHANGELOG.md.

**Cold-start audit (June 2026).** A stratified 25-foil cold-start sweep
(`bench/cold_start_sweep.py`, see `bench/results/REPORT.md`) re-measured these.
On the current branch, **4 of the 5 single-point limitations below now converge
cold** — the entries are retained for history but were already resolved by prior
work. Only **NACA 0008-34 α=−2.6° still fails cold.** The audit also added an
XFOIL-style RMS convergence criterion (Phase 1) and an `acoustic_nan`
failure_mode (a converged aero solve whose noise model returns non-finite OASPL,
common at low Re — formerly a silent blank failure). A residual line search and
ctau equilibrium reseeding were tried and rejected as net-neutral/regressive
(`bench/results/PHASE2.md`, `PHASE3.md`).

- **NACA 0008-34 α=−2.6°** *(still fails cold)*: genuine multi-node BL attractor
  (period-14 cycle). Warm-start from adjacent alphas also fails. No fix found;
  documented as permanent.

- **Cold-start period-2 oscillation** *(documented points converge cold; class
  not eliminated)*: Boeing 737 Midspan α=−3.1°/−3.2°. Verified at the **documented
  condition** (Re=2e6, nCrit=5): both converge cold (35/34 it), as do they at
  nCrit=9. But the transition-sensitivity persists — the neighbour α=−3.0° still
  cold-fails with `transition_front_oscillation`. So the specific documented
  points are no longer stuck, but "others near transition-sensitive points" still
  cycle; `failure_mode="transition_front_oscillation"` signals them.

- **NACA 0012, nCrit=5, α=±2.5°** *(converges cold now, ~6 it)*: was NaN-lock —
  BL Jacobian going singular at iter 3 and staying frozen (`failure_mode="nan_lock"`).

- **NACA 0012, nCrit=5, α=4.7°** *(converges cold now, ~11 it)*: was a multi-node
  BL attractor with ctau-limiter saturation; residual plateaued at ~0.17.

### Forward convergence knob
`solve_coupled` converges on an RMS residual (`resid_rms` in coupled.cpp) against
`param.rtol` (default 1e-6, XFOIL-comparable). `rtol` is plumbed from
`input["rtol"]` (forward-only; AD path unaffected) so it can be swept without
rebuilding. It is now exposed on `OperatingConds(rtol=...)` (default 1e-6); tighten
to ≤1e-10 for AD-vs-FD gradient checks — central FD amplifies converged-state noise
by 1/(2h), so a loose forward solve manufactures spurious FD-vs-AD error (see
`bench/results/FREE_TRANS_VERIFY.md`).
