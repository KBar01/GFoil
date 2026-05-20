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

### Do NOT add Real fields to Param_t<Real> in data_structs_shared.hpp
`Param_t<Real>` is instantiated inside CoDi active tape regions (e.g. inside
`partialOutputspartialInputs` after `tape.setActive()`). Any `Real` field with a
default initializer (`Real foo = X`) will execute a tape assignment statement at
construction time, creating a spurious tape entry and silently corrupting AD
gradients. This was confirmed when adding `ncrithyst` to `Param_t<Real>` shifted
gradient arrays by ~6×10⁻⁷.
Rule: fields needed only by the forward solver belong in the non-template `Param`
struct in `data_structs.h`, not in `Param_t`.

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

## Bug Fixes

### Transition location period-2 limit cycle fix (May 2026) — COMPLETE
Failure mode: at certain operating conditions (confirmed at alpha=−4.9°,
Re=2×10⁶, Ma=0, ncrit=5) the Newton solver entered a stable period-2 limit
cycle. Residual reached ~1.5×10⁻⁴ then oscillated indefinitely between two
states differing only in amplification factor at the last-laminar node.
Transition node index did NOT move — the oscillation was entirely in the amp
values feeding back into the Newton system.

Root cause: in `update_transition.cpp`, the `ilam == ilam0` branch (transition
node unchanged) restored ALL `sa[]` values including the ODE-consistent laminar
amplifications computed by `march_amplification`. This discarded the
march-computed values and reinstated the Newton-perturbed values, which
oscillated with period-2 amplitude ~3.2×10⁻³ at `Is[ilam0]`. The corrupted
amp then fed back into the next Newton step, perpetuating the cycle.

Fix: in the `ilam == ilam0` branch, restore only **turbulent** nodes' ctau from
`sa[]` (undoing march's corruption of ctau at turbulent nodes), but keep the
march-computed laminar amp values which satisfy the eN ODE exactly. One branch
changed in `update_transition.cpp`.

Verified:
- alpha=−4.8° (was already converging): still converges; CL/CD/OASPL differ by
  2–7×10⁻⁶ relative from pre-fix (expected: different convergence path near
  transition changes which amp values feed into the final Newton iterate)
- alpha=−4.9° (previously 2-cycled forever): now converges directly without
  continuation; CL/CD/OASPL match continuation reference to <10⁻⁴ relative

`ncrithyst` parameter: added to `Param` in `data_structs.h` (default 0.2) and
plumbed through the full stack (`run_forward.h`, `run_forward.cpp`, `main.cpp`,
`gfoil_fwd_bindings.cpp`, `OperatingConds` in `inputs.py`, `_build_input_dict`
in `gfoil.py`) but not yet actively used in any computation. Retained as a
future tuning knob for the brief's originally-proposed hysteresis approach.

**CRITICAL CoDi constraint discovered:** do NOT add `Real` fields with default
initializers to `Param_t<Real>` in `data_structs_shared.hpp`. `Param_t<Real>`
is instantiated inside CoDi active tape regions; a new `Real` field creates a
spurious tape entry and corrupts AD gradients silently. `ncrithyst` is therefore
kept only in the non-template `Param` in `data_structs.h`.

Golden files regenerated: prior golden files were stale (GFoil_AD binary
crashed against them). New golden files pass 10/10 at 0.00e+00 error.

Files modified: `src/update_transition.cpp`, `src/include/data_structs.h`,
`src/include/run_forward.h`, `src/run_forward.cpp`, `src/main.cpp`,
`src/gfoil_fwd_bindings.cpp`, `GFoil/inputs.py`, `GFoil/gfoil.py`

---

## Performance Optimisations

### AIC panel geometry precomputation (May 2026) — COMPLETE
Hot path: `build_gamma_codi` in `solve_inv.hpp` runs a 200×200 double loop to
assemble the influence coefficient matrix. Previously each (i,j) iteration
called `panel_info()` in full, recomputing the panel-fixed quantities (t, n, d)
200 times per panel.

Changes made:
- `panel_funcs.hpp`: added `PanelGeom<Real>` struct holding `t[2]`, `n[2]`, `d`
- `panel_funcs.hpp`: added `precompute_panel_geom()` — fills `PanelGeom` from two
  endpoint coordinates; uses `sqrt(dx*dx+dy*dy)` directly
- `panel_funcs.hpp`: added `panel_info_cp()` — fills `PanelInfo` given a
  precomputed `PanelGeom` and a control-point position; copies t/n/d from
  struct, computes x/z/r1/r2/theta1/theta2 only
- `panel_funcs.hpp`: added two new overloads of `panel_linvortex_stream` and
  `panel_constsource_stream` that accept `PanelGeom` instead of endpoint coords
- `solve_inv.hpp`: allocates `PanelGeom<Real> panelGeoms[Ncoords]` on the stack
  before the i-loop; indices 0..Ncoords-2 for body panels, Ncoords-1 for TE
  panel; all three stream calls inside the loop use the precomputed overloads

Deferred: `info.d = norm_t_init` fix in `panel_info()` — NOT applied. The
original `panel_info` is still used in `calc_ue_m.hpp` paths; applying the fix
there would change the CoDi tape accumulation order and shift gradients.
Left as a known minor inefficiency.

Golden files regenerated: `PanelGeom<Real>` stores CoDi active types, so
reordering the tape accumulation causes floating-point associativity shifts
of 3–6×10⁻⁸ relative in gradient arrays. Forward scalars (CL/CD/CM/OASPL)
are bit-for-bit unchanged. New golden files committed.

Files modified: `src/include/panel_funcs.hpp`, `src/include/solve_inv.hpp`

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

### Multiple observer locations (COMPLETE)
### Pybind11 bindings (COMPLETE)
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

### A-weighting (COMPLETE)

### Pybind11 Python bindings
Defer until C++ API is stable (after multiple observers + forced
transition). The restart.json Jacobian state will be passed as a
Python object between run_forward() and run_AD() calls, eliminating
file I/O. See design notes in conversation history.

## Current Work
Multiple observer locations: COMPLETE.
Pybind11 bindings + clean Python API: COMPLETE.
  - `gfoil_cpp.cpython-38-*.so` built into `GFoil/GFoil/` (importable as `from . import gfoil_cpp`)
  - `fwd_run()` returns `FwdResult`; `grad_run(result, ...)` returns `GradResult`
  - No file I/O in pybind11 path; subprocess fallback still available
  - Build: `cmake -B build . -DPYBIND11_PYTHON_VERSION=3.8 -DPYTHON_EXECUTABLE=$(pyenv which python3.8) && cmake --build build --target gfoil_cpp`
Warm-start in pybind11 continuation path: COMPLETE.
A-weighting toggle: COMPLETE.
AIC panel geometry precomputation: COMPLETE.
Transition period-2 limit cycle fix: COMPLETE.
Next task: remaining noise-path optimisations (Estar caching, surf-loop
hoisting), then pyOptSparse integration or forced transition.

### Completed since last CLAUDE.md update
- Transition period-2 limit cycle fix (May 2026):
    - See "## Bug Fixes" section above for full details
    - One branch changed in update_transition.cpp (`ilam == ilam0`: restore
      only turbulent ctau from sa[], keep march-computed laminar amps)
    - ncrithyst plumbed through stack (default 0.2, unused in computation)
    - CRITICAL: ncrithyst NOT in Param_t<Real> — spurious tape entry risk
    - Golden files regenerated (stale; now 10/10 at 0.00e+00)
- AIC panel geometry precomputation (May 2026):
    - See "## Performance Optimisations" section above for full details
    - Golden files regenerated: floating-point associativity shift only
      (3–6×10⁻⁸ relative in gradient arrays); forward scalars unchanged
- A-weighting toggle (May 2026):
    - calc_OASPL gains `const int aWeighting = 0` as final parameter (after WPSjson)
    - Inside the per-observer loop, after TE_noise_outer and before the iObs==0
      cache block, an A-weighting loop multiplies farfieldSpectra[i] by RA²
      using the IEC 61672 formula; all arithmetic stays as Real with
      `static_cast<Real>(precomputed_double)` constants
    - runCode (run_forward.h / run_forward.cpp) gains `int aWeighting = 0`,
      forwarded to calc_OASPL
    - src/main.cpp: both WPSonly branch and normal branch read
      `j.value("aWeighting", 0)` and pass through
    - gfoil_fwd_bindings.cpp: reads aWeighting from input dict
      (inp.contains guard for backward compat), passes to runCode
    - gfoil_ad_bindings.cpp: same, passes to partialOutputspartialInputs
    - srcAD/include/ADfuncs.hpp: partialOutputspartialInputs gains
      `int aWeighting = 0`, forwarded to calc_OASPL
    - GFoil/inputs.py: Acoustics gets `aWeighting: bool = False`
    - GFoil/gfoil.py: _build_input_dict includes `"aWeighting": int(acoustics.aWeighting)`
    - AD differentiates through A-weighting correctly (pure Real arithmetic,
      no tape detachment); golden files unchanged (aWeighting=0 by default)
- Warm-start fix in pybind11 path (May 2026):
    - runCode gains `const RestartState* warmStart = nullptr` parameter
    - When non-null, initialises glob.U and vsol.turb from in-memory state
      instead of reading restart.json
    - run_forward_py gains optional `prev_jacobian` argument
    - _call_forward in gfoil.py passes prev_result states/turb through
    - standard_run tracks last_converged and passes it as warm start
      in the forward-stepping loop
    - Subprocess fallback: writes restart.json from prev_result if available
    - Binary path unchanged: fromRestart=1 still reads restart.json
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
  A-weighting feature:   +~60 lines net (new param plumbing + IEC loop)
  ─────────────────────────────────────
  Total net reduction:  ~-1219 lines
  Duplicate function definitions eliminated: 13
  Dead code removed: Trans struct + 5 functions + forced-trans plumbing