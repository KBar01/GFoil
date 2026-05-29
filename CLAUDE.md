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
in `gfoil.py`). Activated as a single-node retreat gate in `update_transition.cpp`
(see "Bug Fixes — ncrithyst" entry). The regression test uses `ncrithyst=0`.

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

### Ctau freeze cycle-detection mechanism (May 2026) — COMPLETE
Failure mode: at transition-sensitive operating points, the Newton solver
can enter a period-N limit cycle after transition stabilises. Residual
oscillates indefinitely above tolerance even though the transition node
index (ilam) is no longer moving.

Mechanism: when (a) ilam has been stable for ≥8 consecutive iterations
AND (b) the L2 residual has not improved by 2× over the last 8 iterations
(circular residual buffer in `coupled.cpp`), a per-surface `ctau_freeze`
flag is set. While frozen, `update_transition.cpp` averages the ctau at
the first turbulent node (Is[ilam0+1]) between its current Newton-updated
value and the previous Newton-updated value, storing the pre-averaging value
so a period-2 cycle collapses after exactly 2 freeze iterations.

Co-activation: once either surface freezes, the partner surface is also
frozen if it has had ≥4 consecutive stable-ilam iterations (prevents
cross-coupling oscillations when one surface cycles while the other is frozen).

Release: freeze is maintained until ilam moves (transition shifts) or the
solver converges. Releasing early resets the averaged history and restarts
the same cycle.

CoDi compatibility: freeze detection uses only plain doubles (resid_buf,
prev_ctau, ilam counters) — no getValue() on design variables. The single
getValue() call inside the freeze block reads ctau at Is[ilam0+1] and writes
it back as Real((prev+curr)*0.5), which detaches that node's ctau from the
CoDi tape only when freeze is active. For the regression golden case (alpha=2°,
ncrit=5) the freeze never activates (converges in 10 iterations), so the AD
computation is unaffected.

Effectiveness: the freeze successfully resolves single-node ctau oscillations
(period-2 cases such as Boeing 737 Midspan α=−3.1°/−3.2°) by collapsing the
cycle after exactly 2 averaging iterations. It cannot break multi-node coupled
BL attractors of the type originally diagnosed at NACA 0008-34 α=−2.6° (see
Known Limitations); that case was subsequently resolved by the ncrithyst
hysteresis activation (see ncrithyst entry below).

Files modified: `src/coupled.cpp`, `src/update_transition.cpp`,
`src/include/main_func.h`

Regression test: golden files regenerated after restoring input.json to
alpha=2°, ncrit=5.0, NACA 0012 (n0012_sharp.dat), observer [0, 3, 0.5].
10/10 at 0.000e+00. Note: OASPL in this test depends only on observer x and z
(not y) since newAmiet.hpp uses mid-span formulation (y ignored in S0).

### Transition-node jump cap (May 2026) — COMPLETE
Failure mode: cold-start Newton solve at certain alpha/nCrit combinations
(confirmed at alpha=7.3°, Re=2×10⁶, nCrit=6) exhausted all 60 outer
iterations without converging. Residual plateaued at ~1.1–1.3, omega locked
at 0.001–0.02 for 57 of 60 iterations.

Root cause: at alpha=7.3° the cold-start `init_BL` march left the amplification
factor 2.34 below ncrit at the transition node (vs 1.44 at alpha=6.92° which
converges). The large amp residual drove the first few Newton steps to overshoot
the transition location from node 33 to node 63 in 4 iterations (30-node jump).
This filled 30 nodes with inconsistent BL state (laminar theta/delta_star,
turbulent ctau from `get_cttr`), driving the ctau limiter in `update_state`
(`max ctau change = 0.05/step`) to omega~0.001 for the remainder of the run.

Fix: in `update_transition.cpp`, immediately after `march_amplification` returns
`ilam`, cap the forward advance at 3 nodes per Newton iteration:
```cpp
constexpr int max_transition_jump = 3;
if (ilam < ilam0) {
    ilam = std::max(ilam, ilam0 - max_transition_jump);
}
```
The cap only applies to the `ilam < ilam0` branch (transition moving forward /
more surface becoming turbulent). The `ilam > ilam0` branch (retreat) is
unchanged — retreating transition removes turbulent nodes, which is benign for
the ctau limiter.

Verified:
- alpha=7.3°, nCrit=6 cold-start: converges in 1 call, ~0.58s (was 7 calls,
  3.45s with backstepping)
- alpha=6.5°, nCrit=6 cold-start: converges in 1 call, ~0.55s (was failing cold,
  needed backstepping: 4 calls, 1.46s)
- alpha=6.92°, 8.0°, 9.0°: unchanged
- alpha=−4.9° (period-2 fix case): unchanged
- Regression golden (alpha=2°, nCrit=5): bit-for-bit identical (cap never
  triggers when transition is stable)

Instrumentation added (guarded by `GFOIL_DEBUG=1` env var, silent by default):
- `src/coupled.cpp`: per-iteration print of L2 residual, omega, and last-laminar
  node index for both surfaces, plus INIT state print from `init_BL.cpp`
- `src/update_state.cpp`: return type changed `void → Real` (returns omega) so
  `coupled.cpp` can print it without duplicating limiting logic
- `src/include/main_func.h`: declaration updated to match

Files modified: `src/update_transition.cpp`, `src/coupled.cpp`,
`src/update_state.cpp`, `src/include/main_func.h`, `src/init_BL.cpp`

### ncrithyst transition hysteresis activation (May 2026) — COMPLETE
`ncrithyst` (default 0.2) was plumbed through the full stack when the period-2
fix was implemented but left inactive (unused in any computation).

Design constraint (derived from testing): applying hysteresis as a threshold
offset inside `march_amplification` shifts the ODE solution — physically wrong.
Applying it as a gate on BOTH advance and retreat in `update_transition` breaks
convergence because `amp_break` is always only slightly above `ncrit` (set by
the ODE), so the `amp_break > ncrit+ncrithyst` advance gate fires at every
normal convergence step.

Correct implementation: `march_amplification` uses `param.ncrit` exactly
(unchanged). The gate is applied only to **single-node retreats** in
`update_transition` — the minimum change that damps spurious 1-node
laminar-recovery oscillations without affecting any advance direction:

```cpp
// march_amplification: unchanged ncrit threshold; amp_break output added
if (U2[2] > param.ncrit) {
    if (amp_break) *amp_break = U2[2];
    break;
}

// update_transition: gate only 1-node retreats
} else if (ilam > ilam0 && (ilam - ilam0 == 1) && (ilam0 + 1 < nSurfPoints)) {
    Real amp_first_turb = glob.U[colMajorIndex(2, Is[ilam0+1], 4)];
    if (amp_first_turb >= param.ncrit - param.ncrithyst) {
        ilam = ilam0;  // hysteresis: suppress spurious 1-node retreat
    }
}
```

Rationale: advance is ODE-authoritative (march already returns the physically
correct transition node). Only retreat needs damping: a single-node laminar
recovery oscillation where amp barely drops below ncrit is suppressed when the
node's virtual amp (from march) is still within the hysteresis band. Multi-node
retreats needed for convergence from an off-equilibrium initial state are always
unrestricted.

With `ncrithyst=0`: gate condition `amp_first_turb >= ncrit - 0 = ncrit` is
essentially never met (amp < ncrit by march construction), giving bit-for-bit
identical results to pre-ncrithyst. The regression test uses `ncrithyst=0` in
`input.json` for this reason.

Verified (NACA 0012 unless noted, `ncrithyst=0` unless noted):
- alpha=2°, ncrit=5, ncrithyst=0 (regression golden): bit-for-bit identical to
  pre-ncrithyst. 10/10 at 0.000e+00 against restored golden files.
- alpha=2°, ncrit=5, ncrithyst=0.2: converges; CL shifts by ~1×10⁻³ rel (retreat
  gate alters convergence path slightly near ncrit=5). Expected physical effect.
- alpha=−4.9°, ncrit=5, ncrithyst=0.2: converges; CL matches ncrithyst=0 to ~2×10⁻¹¹
  rel (retreat gate never fires for this case).
- alpha=7.3°, ncrit=6, ncrithyst=0.2: converges; CL shifts by ~3×10⁻³ rel (retreat
  gate fires on some convergence iterations near this ncrit).
- NACA 0008-34, alpha=−2.6°, ncrit=9: fails both with and without ncrithyst.
  This case is the genuine multi-node coupled BL attractor (period-14 cycle).
  The correct hysteresis implementation does not resolve it (see Known Limitations).

Full 12-aerofoil sweep (Re=2e6, Ma=0, ncrit=9, ncrithyst=0.2, 1452 cases):
  Cold converged: 83.8%  (1217/1452)
  Total converged: 98.8% (1434/1452)
  Failures: 18

Files modified: `src/update_transition.cpp`, `tests/golden/*.json`

---

## Known Limitations

### Cold-start period-2 oscillation (May 2026) — ACCEPTED LIMITATION

Certain alphas near transition-sensitive operating points fail to converge
on a cold start but converge correctly via warm-start continuation from a
nearby alpha. Diagnosed cases:
  Boeing 737 Midspan: α = −3.2°, −3.1°
  NACA 2414:          α = −5.1° (pre-existing pybind11 crash, separate issue)

Root cause: a large transition retreat at iter 3 (up to 78 nodes) leaves
ctau values at the transition front far from equilibrium. The 0.05/step
ctau limiter in update_state prevents convergence within 60 iterations.
A period-2 oscillation develops once transition stabilises at the correct
location: the Newton Jacobian makes a sign error in the ctau correction
at the transition-front node, creating alternating over/undershoots.

Approaches investigated and rejected:
- Laminar amp cap in init_BL.cpp: correct invariant, kept, but does not
  fire during the final oscillation phase (amp stays below ncrit then)
- Symmetric retreat cap in update_transition.cpp: caused regression on
  Boeing 737 Outboard α=−0.7° (24-node retreat physically necessary there)
- Residual-adaptive ctau limit in update_state.cpp: relaxing the limit
  near convergence made oscillation worse (larger ctau step → larger
  overshoot the next iteration); reverted

Mitigation: `failure_mode = "transition_front_oscillation"` is returned on
FwdResult so callers can detect the pattern and retry with a warm start.
The sweep script's continuation logic already handles these cases correctly.

### Ctau multi-node coupled BL attractor — NACA 0008-34 α=−2.6° (May 2026) — GENUINE FAILURE

**Distinct from cold-start period-2 oscillation** (see above). The cold-start
oscillation cases (Boeing 737 Midspan α=−3.1°/−3.2°) are characterised by:
  - A large initial transition retreat that leaves ctau far from equilibrium
  - Single-node ctau cycling at the transition front after ilam stabilises
  - Resolution by warm-start continuation from a neighbouring alpha

NACA 0008-34 α=−2.6° is categorically different:
  - Warm starts from α=−2.7° AND α=−2.5° both fail identically — this is
    confirmed NOT cold-start sensitivity
  - The transition node settles normally (no large initial retreat), then
    a stable period-14 limit cycle develops within the final Newton approach phase
  - The cycle is a multi-node coupled BL attractor: multiple turbulent nodes near
    the transition front oscillate simultaneously through the Newton Jacobian, not
    a single-node ctau oscillation at Is[ilam0+1]
  - The ctau freeze activates correctly on both surfaces and reduces cycle
    amplitude from ~1×10⁻² to ~5×10⁻⁴, but the minimum residual ~5×10⁻⁴ does
    NOT decrease across freeze iterations — it is a true stable attractor
  - Single-node ctau intervention at Is[ilam0+1] is insufficient because the
    coupling involves BL state at multiple coupled nodes simultaneously
  - The ncrithyst retreat hysteresis does not resolve it: the period-14 cycle
    involves a multi-node coupled attractor, not a single-node retreat oscillation

Root cause: at this alpha/ncrit combination the Newton Jacobian simultaneously
makes a sign error in the ctau correction at multiple turbulent nodes near the
transition front, creating correlated alternating overshoots that single-node
averaging cannot decouple.

Approaches investigated and rejected:
- get_cttr anchor for transition-front ctau: oscillates because get_cttr reads
  the already-oscillating BL state (theta/ds/ue) from glob.U
- ctau averaging over 3 nodes (instead of 1): made things worse (minima ~4×10⁻³
  vs ~5×10⁻⁴ for 1-node averaging)
- 1e-3 residual release of freeze: re-triggers the cycle immediately
- ncrithyst hysteresis activation: correctly gates single-node retreats but
  cannot decouple the multi-node Jacobian coupling

`failure_mode = "transition_front_oscillation"` is returned so callers can
detect it, but warm-start continuation cannot rescue this case.

### pybind11 in-process segfault fix (May 2026) — COMPLETE

**Symptom**: `convergence_sweep.py` (using the pybind11 in-process path) crashed
deterministically around case 77 with a SIGSEGV inside
`Eigen::SparseLU::solve()` / `MappedSuperNodalMatrix::solveInPlace()`.
The subprocess path (standalone `GFoil_fwd_codi`) did not crash because each
call is an isolated process.

**Root cause**: In `solve_sys_sparse` (`src/include/sparselinsolve.hpp`), a call to
`lu.compute(A)` with a matrix containing NaN/Inf entries fails silently with
`info = NumericalIssue`. The subsequent `lu.solve(b)` then dereferences a null
`supToCol()` pointer from the uninitialised factorisation internals, causing
the SEGV.

**Origin of NaN in the Jacobian**: During Newton iterations for certain
aerofoil/alpha combinations, BL residual station computations produce NaN in
the Jacobian blocks at specific nodes (typically consistent across cases: always
the same consecutive-node pair). The NaN is transient — the Newton solver
recovers over subsequent iterations via `stagpoint_move`/`update_transition`
state updates even when `dU=0` — and does not indicate a fundamental physics
failure. It is a pre-existing BL-solver behaviour, not caused by these changes.

**Fix** (`src/include/sparselinsolve.hpp`):
- Added `#include <cstdio>` and `#include <cmath>`
- After `lu.compute(A)`, check `lu.info() != Eigen::Success`
- On failure: set `glob.dU[i] = 0` for all i (no Newton step this iteration)
  and return early; tape registration path is skipped correctly
- Diagnostic (NaN count, location, out-of-bounds indices) printed only when
  `GFOIL_DEBUG=1`, silent in production

**Effect**: The Newton loop skips dU application for NaN iterations; the BL
solver recovers through transition/stagnation point state updates and converges
normally on subsequent iterations. No loss of convergence rate in practice:
sweep statistics improved to 84.6% cold / 98.8% total (≥ prior 83.8%/98.8%).

Regression: 10/10 at 0.000e+00 (golden test case never triggers NaN).

File modified: `src/include/sparselinsolve.hpp`

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
Test input:   tests/input.json  (committed; NACA 0012, alpha=2°, ncrit=5, ncrithyst=0)
              regression_test.py copies this to repo root before running binaries.
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
Transition-node jump cap: COMPLETE.
Ctau freeze cycle-detection mechanism: COMPLETE (see Bug Fixes).
  - Detects and partially mitigates period-N ctau oscillations at stable
    transition fronts.
ncrithyst hysteresis activation: COMPLETE (see Bug Fixes).
  - march_amplification uses ncrit exactly; only single-node retreats gated.
  - Regression test uses ncrithyst=0 (bit-for-bit identical to pre-ncrithyst).
  - Sweep (ncrithyst=0.2, ncrit=9): 83.8% cold / 98.8% total / 18 failures.
  - NACA 0008-34 α=−2.6° remains a genuine failure (multi-node attractor).
pybind11 in-process segfault fix: COMPLETE.
  - SparseLU::compute NaN failure now handled gracefully (dU=0 skip iteration)
  - Diagnostic available under GFOIL_DEBUG=1; silent by default
  - Sweep: 84.6% cold / 98.8% total; regression: 10/10 at 0.000e+00
Next task: pyOptSparse integration or forced transition (cold-start
  oscillation for 2+ specific alphas is documented as accepted limitation).

### Completed since last CLAUDE.md update
- Transition period-2 limit cycle fix (May 2026):
    - See "## Bug Fixes" section above for full details
    - One branch changed in update_transition.cpp (`ilam == ilam0`: restore
      only turbulent ctau from sa[], keep march-computed laminar amps)
    - ncrithyst plumbed through stack (default 0.2, activated in subsequent change)
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
- Transition-node jump cap (May 2026):
    - See "## Bug Fixes" section above for full details
    - 3-node-per-iteration cap on forward transition advance in
      update_transition.cpp (`ilam < ilam0` branch only)
    - update_state return type changed void→Real to expose omega for
      GFOIL_DEBUG=1 instrumentation; main_func.h declaration updated
    - GFOIL_DEBUG=1 env var enables per-iteration residual/omega/ilam
      diagnostics from coupled.cpp and init state from init_BL.cpp
    - Golden files unchanged (cap never triggers at golden alpha=2°)
- Ctau freeze cycle-detection mechanism (May 2026):
    - See "## Bug Fixes" section above for full details
    - Circular residual buffer (8 entries) + stable_ilam counters in coupled.cpp
    - freeze flag activates when ilam stable ≥8 iters + residual not improving 2×
    - In ilam==ilam0 branch: when frozen, average ctau at Is[ilam0+1] between
      current and previous Newton values (store pre-averaging value)
    - Co-activation: if one surface frozen + partner stable ≥4 iters, freeze both
    - Reduces cycle amplitude; NACA 0008-34 α=−2.6° attractor remains genuine failure
    - CoDi safe: getValue() only at Is[ilam0+1] when freeze is active
    - Golden files regenerated: input.json restored to alpha=2°, ncrit=5,
      NACA 0012 (n0012_sharp.dat), observer [0,3,0.5]; 10/10 at 0.000e+00
    - Files modified: src/coupled.cpp, src/update_transition.cpp,
      src/include/main_func.h, tests/golden/*.json
- ncrithyst hysteresis activation (May 2026):
    - See "## Bug Fixes — ncrithyst" entry above for full design rationale
    - march_amplification: ncrit threshold unchanged; Real* amp_break output added
    - update_transition: single-node retreat gate added (ilam > ilam0 and
      ilam-ilam0==1): blocks retreat if amp_first_turb >= ncrit-ncrithyst
    - Advance direction: no gate — march is ODE-authoritative for advance
    - Regression: input.json uses ncrithyst=0 for bit-for-bit identity; 10/10
    - Sweep (ncrithyst=0.2, ncrit=9): 83.8% cold / 98.8% total / 18 failures
    - NACA 0008-34 α=−2.6° unchanged: remains genuine failure
    - Files modified: src/update_transition.cpp, tests/golden/*.json
- pybind11 in-process segfault fix (May 2026):
    - Root cause: NaN/Inf BL Jacobian entries → SparseLU::compute fails → null
      supToCol() → SEGV in lu.solve(b). Pre-existing NaN in BL residual blocks
      at certain operating conditions; previously masked in subprocess path.
    - Fix: check lu.info() after lu.compute(A); on failure set dU=0, return early
    - Diagnostic (NaN count, location) gated on GFOIL_DEBUG=1; silent by default
    - Sweep improved: 84.6% cold / 98.8% total; regression: 10/10 at 0.000e+00
    - Files modified: src/include/sparselinsolve.hpp

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