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

**Never include both in the same TU** — `norm2` redefinition results.
(This was the error seen in the May 2026 build.)

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
When adding support for multiple observer locations, pass them as:
  const Real* obsX, const Real* obsY, const Real* obsZ
or as:
  const Real (&obsX)[N]
Never as double* or std::vector<double>.

---

## Shared solver template headers

`src/include/solver_funcs.hpp` — included by both build targets.
Contains duck-typed template implementations; each fwd *.cpp delegates
via a one-line non-template wrapper. main_func.hpp includes it directly.

`src/include/panel_funcs.hpp` — already shared; `inviscid_velocity` and
`dvelocity_dgamma` added as `template<Real,FoilT>` in May 2026.

## Dead / Stale Files

- `src/noise_includes/amiet.h` — deleted (May 2026); was not included anywhere.

---

## Audit History (May 2026)

**newAmiet.hpp** (921 → 402 lines):
- Audited against R&M (2005) mid-span formulation
- Removed subcritical branch and all unused helpers
- Removed params from `TE_noise_outer`: rho, nu, Ky
- Removed params from `Radiation_integral_total`: Ky, K_2_bar, beta
- `Fresnel_int()` removed; E(x) now computed as conj(E*(x))

**Build status (May 2026):**
- `GFoil_fwd_codi`: builds and links cleanly ✓
- `GFoil_AD`: builds and links cleanly ✓
- All newAmiet.hpp CoDi issues (fabs→abs, hypot, using-declarations) resolved.

---

## Testing / Regression

Golden-output regression test covers CL, CD, CM, OASPL (forward) and
all three gradient arrays + alpha scalars (AD).

  # First-time setup (already done — golden files committed)
  python3 tests/regression_test.py --create-golden

  # Run after every refactoring step
  python3 tests/regression_test.py --test

  # Rebuild then test
  python3 tests/regression_test.py --build --test

Golden files: tests/golden/fwd_scalars.json, ad_scalars.json,
              ad_gradients.json
Tolerance: 1e-8 relative.

Note: both binaries exit with code 1 even on success. The test script
detects success by checking that the output file was written/updated,
not by exit code. Do not change this logic without also fixing the
binary exit codes.

## Deduplication progress

### Completed (Easy tier — all in solver_funcs.hpp or panel_funcs.hpp)
- init_thermo            → solver_funcs.hpp  template<OperT,ParamT,GeomT>
- space_wake_nodes       → solver_funcs.hpp  template<Real,FoilT,WakeT>
- identify_surfaces      → solver_funcs.hpp  template<IsolT,VsolT>
- set_wake_gap           → solver_funcs.hpp  template<FoilT,IsolT,VsolT>
- rebuild_ue_m           → solver_funcs.hpp  template<FoilT,WakeT,IsolT,VsolT>
- stagpoint_find_impl    → solver_funcs.hpp  template<bool,GammaT,VarT,FoilT,WakeT>
- calc_force             → solver_funcs.hpp  template<OperT,GeomT,ParamT,FoilT,GlobT,PostT>
- build_wake_impl        → solver_funcs.hpp  template<FoilT,GeomT,OperT,IsolcT,WakeT>
- inviscid_velocity      → panel_funcs.hpp   template<Real,FoilT>
- dvelocity_dgamma       → panel_funcs.hpp   template<Real,FoilT>

### Completed (Medium tier)
- stagnation_state kernel  → solver_funcs.hpp  stagnation_state_impl<Real>
- stagpoint_move core      → solver_funcs.hpp  stagpoint_move_impl<...>
- Ue residual kernel       → solver_funcs.hpp  ue_residual_kernel<...>
- 9 structs (AD side)      → data_structs_shared.hpp  template<typename Real>
                             AD data_structs.hpp now 54 lines (was 241)

### Overall deduplication results
  Lines removed:                     1473
  Lines added (shared headers):      1057
  Net reduction:                     -416 lines
  Duplicate function defs eliminated: 13

### Known deferred items
- src/include/data_structs.h: fwd side still defines its own 9 structs
  because struct Foo; forward declarations in fwd headers are incompatible
  with template type aliases. Requires a broader include-order refactor.
- build_glob_RV / build_glob_RV_AD: Hard — Jacobian fill interleaved at
  every station; only ~35% logic shared. Do not attempt without careful
  planning.

## Planned Features

### Multiple observer locations (next to implement)
Average OASPL across N observers using power averaging:
  OASPL_avg = 10 * log10( (1/N) * sum_i( 10^(OASPL_i / 10) ) )

The averaged OASPL is the single scalar that the AD differentiates,
so the gradient pipeline (dOASPL/dy, dOASPL/dalpha) is unchanged.

Design:
- calc_OASPL signature changes from scalar (X,Y,Z) to arrays
  (obsX[], obsY[], obsZ[], nObs) — all typed as Real
- WPS computation (calc_WPS) is observer-independent — compute ONCE
  before the observer loop (it only depends on BL states)
- TE_noise_outer is called once per observer inside the loop
- Frequency integration and dB conversion done per observer
- Power average over all observers at the end
- JSON reading in main.cpp: if X/Y/Z are JSON arrays use them
  directly; if scalars wrap in single-element array — backward
  compatible with existing input.json
- inputs.py Acoustics dataclass: accept observerXYZ as shape (3,)
  for single observer or (N,3) for N observers; normalise to (N,3)
  internally; pass as JSON arrays always

Key constraint: ALL observer coordinate arrays must be typed as
Real throughout — never cast to double or pass as double*.

### Forced transition (future branch)
Re-implement cleanly after multiple observer support is complete.
Infrastructure was removed in commit 1065ed2.

### A-weighting
Post-process farfieldSpectra with A-weighting curve before
integration. Apply inside calc_OASPL after TE_noise_outer call,
before the frequency integration loop.

### Pybind11 Python bindings
Defer until C++ API is stable (after multiple observers + forced
transition). The restart.json Jacobian state will be passed as a
Python object between run_forward() and run_AD() calls, eliminating
file I/O. See design notes in conversation history.

## Current Work
Next task: implement multiple observer locations.
See "## Planned Features" above for the full design spec.
Key constraint to remember: ALL types must be Real — no double* for
observer coordinates. See "## CoDi Type Rules (EXTENDED)".

### Completed since last CLAUDE.md update
- newAmiet.hpp physics fixes (supervisor review):
    - Fresnel_int_conj renamed to Estar throughout
    - G_e coefficient corrected: sqrt(0.5*k/D) not sqrt(k/D)
    - Mid-span S0: sqrt(x^2 + beta^2*z^2), y term removed
    - b_half = b (semi-chord passed directly, not c/2.0)
    - std::fabs → std::abs throughout (CoDi compatibility)
    - std::hypot → std::sqrt(a*a + b*b) (CoDi compatibility)
    - using std::complex / using std::exp removed from file scope
    - TODO comments added near denominator clipping in G_c, G_d, G_e
    - R&M equation references added to all function comment blocks
- Note: golden files were regenerated after these changes (physics
  change — OASPL values shifted). Commit includes new golden files.

### Running totals (all refactoring to date)
  Easy tier dedup:       -307 lines net
  Medium tier dedup:     -109 lines net
  Forced trans removal:  -813 lines net
  newAmiet.hpp cleanup:  ~-50 lines net (fabs, using decls, comments)
  ─────────────────────────────────────
  Total net reduction:  ~-1279 lines
  Duplicate function definitions eliminated: 13
  Dead code removed: Trans struct + 5 functions + forced-trans plumbing