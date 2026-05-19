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

### Deferred (Medium difficulty)
- stagnation_state       fwd adds Ust_U[32]+Ust_x[8] Jacobian outputs; AD omits
- stagpoint_move         different preambles (fwd scans glob.U; AD receives currStag)
- Ue residual kernel     fwd adds Jacobian fill + sparse solve on top of shared kernel
- struct unification     Geom,TE,Foil,Wake,Oper,Post,Param,Trans,Vsol all field-identical
                         → template<typename Real> structs in data_structs_shared.hpp

### Deferred (Hard — do not attempt)
- build_glob_RV          Jacobian fill interleaved at every station; ~35% shared logic

## Current Work
About to tackle Medium difficulty deduplication items.
Run python3 tests/regression_test.py --test after every step.