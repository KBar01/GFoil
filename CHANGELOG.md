# GFoil Changelog

Detailed narrative of all significant changes, bug fixes, and investigations.
Standing rules and architecture belong in CLAUDE.md; this file is the
chronological record.

---

## Warm-restart stagnation reindex fixed (June 2026)

**What.** Fixes the second defect scoped in the warm-restart entry. On a warm
entry `runCode` (`run_forward.cpp`) ran `identify_surfaces`/`set_wake_gap`/
`calc_ue_m` from the INVISCID `stagpoint_find`, then a single viscous
`stagpoint_move` whose sign-scan is seeded from `isol.stagIndex`. Left at the
inviscid value it landed one node off the donor's converged stag (e.g. [81,82] vs
donor [82,83]), reindexing the BL stations: a same-alpha restart that should
accept in ~0-2 iterations took 8 (entry BL_rms ~0.6). The donor's converged
`stagIndex` was stored in `RestartState.stag` but unused at load.

**Fix.** `RestartState.stag` is now plumbed through the pybind warm path
(`_call_forward` adds `stag` to `jac_in`; `gfoil_fwd_bindings.cpp` reads it) and,
on warm entry, seeds `isol.stagIndex` before `stagpoint_move`. The sign-scan then
starts from the donor bracket, finds the loaded states consistent with it (no
move), and `stagpoint_move`'s `identify_surfaces` rebuilds `Is`/`distFromStag` to
the donor configuration. `wgap` is already donor-identical (keyed to the inviscid
stag, identical for the same geometry+alpha — confirmed in the Part B
instrumentation). `GFOIL_NOSTAGSEED` keeps the pre-fix behaviour for A/B.

**Mechanism (GFOIL_DEBUG, `bench/stagseed_probe.py`).** EPPLER 399 nCrit=5, warm
6.5<-6.5 (donor converges at stag [82,83]): NOSTAGSEED → post-move [81,82], entry
BL_rms 0.621, Is 82/118, it=8; seed → post-move [82,83], entry BL_rms 7.7e-09, Is
83/117, xi 0.992173 (all matching the donor's converged values), it=1. The fix
reproduces the donor configuration rather than inventing one. `stagpoint_move`
rebuilds `Is` every iteration, so the donor's converged Is tracks its converged
stag; the seed reproduces it. For warm 5.0<-5.5 the inviscid stag already equals
the donor [83,84], so the seed is a no-op there.

**Gates.** (1) Regression **22/22 bit-identical** (golden is cold; warm seed path
not taken). (2) warm 6.5<-6.5 **8 → 1** iteration; warm 5.0<-5.5 stays a real
re-solve under Fix B (it=26, nonzero) with cold-truth CL 1.27228 — the
stale-accept fix remains effective. (3) 133-case rescue (`bench/rescue_cost.py`,
seed vs NOSTAGSEED, same binary): rescue rate **94/133 → 94/133** (unchanged);
converged-iteration cost over the 94 commonly-rescued cases **5157 → 4408
(−14.5 %)**, 70 cheaper / 8 marginally costlier (largest: AH 79-100 C +8/nc9
142→87, GOE 346 +8/nc5 59→19). (4) 840 cold sweep vs the post-Part-C baseline:
**0 conv / 0 failure_mode / 0 newton_iteration diffs** (warm path dead on cold
runs). AD path untouched.

---

## Warm-restart stale accept: driver guard + solver root cause + scatter-add (June 2026)

Three sequential changes; full evidence in `bench/results/REPORT.md`
(addenda Part A/B/C). The forward solver and the Python driver mishandled warm
(restart) entries; Part C is an unrelated, bit-identical performance change to the
Newton linear solve. The AD path (`srcAD/`, `gfoil_ad_bindings.cpp`,
`ADfuncs.hpp`) was not touched.

### Part A — Python guard against the warm-restart stale accept (`GFoil/gfoil.py`)

**What.** `run_forward(inp, restart)` could return `converged` with
`newton_iterations == 0` and the DONOR's solution unchanged after an alpha change
(EPPLER 399, nCrit=5: warm 5.0←5.5 returned the 5.5 donor CL 1.3206 as
"alpha=5.0", +3.7% error). `standard_run`'s forward-stepping loop now rejects an
it==0 accept at a changed alpha (`_is_stale_warm_accept`: converged AND it==0 AND
|Δalpha|>1e-12 — geometry/Re/nCrit are invariant within `standard_run`, so alpha
is the only varying input) and retries the same alpha cold. A genuine same-alpha
0-iteration restart is left allowed.

**Measurements.** `bench/eppler399_guard.py`: warm 5.0←5.5 now returns the
cold-truth CL 1.27228 (|Δ|=2.8e-6) via a fresh cold solve; same-alpha 6.5←6.5 not
rejected. `bench/rescue_guard.py` re-ran the 133 cold-failure rescues with the
guard: corrected rescue count **97 → 87** (13 prior "rescues" were stale accepts;
3 newly rescued via healthier continuation; guard fired in 32 runs). CSV gains a
`rescue_conv_guarded` column (`rescue_conv` kept). Golden regression bit-identical
(no C++, golden uses cold `run_forward`).

### Part B — solve_coupled entry-accept root cause (`src/coupled.cpp`, `run_forward.cpp`, `main_func.h`)

**Mechanism (confirmed, GFOIL_DEBUG, `bench/b1_instrument.py`).** `solve_coupled`
tests `resid_rms(glob.R, Rsize)` with `Rsize = 3*(Ncoords+Nwake)` — the BL-station
rows only. The ue-coupling rows (`3*Nsys..4*Nsys-1`), where alpha enters via
`ue_residual_kernel` (called from `solve_glob` AFTER the test), never participate.
Warm 5.0←5.5 ENTRY: BL_rms=**1.8e-07** (< rtol 1e-6) but ue_rms=**0.042**
(>> rtol) → instant it=0 accept. Cold 5.0 control: both large.

**Same-alpha it=8 anomaly — a SECOND, distinct defect (scoped, not fixed).** Donor
cold-6.5 *converges* at `stagIndex=[82,83]`; warm 6.5←6.5 *enters* at `[81,82]`
(Is sizes/`distFromStag` shifted, entry BL_rms=0.62). On a warm restart `runCode`
runs the inviscid `stagpoint_find` then a single viscous `stagpoint_move`, which
lands one node off the donor's converged stag, reindexing the BL stations — so the
BL rows are no longer the donor's converged residuals. `RestartState.stag` is
stored but never re-imposed on warm entry. This is why the original 0.5-deg vs
same-alpha asymmetry exists (instant accept only when the re-entry stag happens to
match). `wgap` was identical (same geometry+alpha → same inviscid stag). A fix
would re-impose the donor stag (+`set_wake_gap`) on warm entry.

**Fixes (both retained behind `GFOIL_BFIX`; default Fix B).** Fix A
(`GFOIL_BFIX=A`): honest criterion `max(BL_rms, ue_rms) < rtol` (a max, not a
pooled 4*Nsys RMS, to preserve the BL rows' per-equation rtol). Fix B (default):
skip the iteration-0 accept on a warm entry (`warmEntry` threaded from `runCode`),
forcing ≥1 Newton iteration; cold paths bit-identical by construction.
`GFOIL_BFIX=off` reproduces the defect.

**Measurements.** EPPLER warm 5.0←5.5: off it=0/1.3206 (stale); A it=26/1.27228;
B it=26/1.27228 (`bench/b3_eppler.py`). Golden regression: Fix B 22/22
bit-identical; Fix A also 22/22 bit-identical (on well-behaved solves the ue rows
are already converged when the BL rows are). 840-case cold sweep
(`bench/cold_only_sweep.py`) vs baseline: BOTH fixes 0 conv / 0 failure_mode /
0 newton_iteration diffs, 707/840, median 16/p90 41/max 59, none over the 60 cap —
so no newly-failing case under Fix A and the B.3(iii) spot-check has no candidates
(the old criterion never declared convergence on an unconverged ue system on this
cold grid). Wall-time (golden, 30 reps): off≈90 ms, B≈90 ms (no-op on cold),
A≈92 ms (~2%, the per-iteration ue-kernel).

**Decision.** Fix B active by default (sufficient, provably bit-identical on cold,
zero-cost); Fix A retained behind the switch (the more principled criterion,
empirically bit-identical here but not universally guaranteed). Baseline contract
unchanged.

### Part C — scatter-add replaces per-iteration setFromTriplets (`src/include/sparselinsolve.hpp`)

**What.** `setFromTriplets` + sparse copy was ~17% of forward time. Keyed to the
existing FNV pattern hash: on a pattern change, build A with `setFromTriplets`,
`analyzePattern`, and record `slot[k]` (the `valuePtr()` index for triplet k, via
binary search of the compressed CSC) and `is_first[k]`. On an unchanged pattern,
skip `setFromTriplets`: the first triplet per slot ASSIGNS (seeds verbatim),
duplicates ADD in ascending-k order — the cached matrix is reused for factorize,
solve, and the unchanged external-function adjoint (`data->A`).

**Bit-identity (mandatory, GFOIL_CVERIFY).** Initial mismatch was signed zero:
the triplet list is duplicate-free here (`nonZeros()==nnz`; `addColumnValues`
pre-sums in place), so `setFromTriplets` stores `-0.0` verbatim while
`fill(0.0)+=` yields `+0.0`. Fixed by assign-first. After: **0 memcmp failures
across 297 solves** (golden + 5 slow cases, incl. golden AD recording);
`dup_iters=0` everywhere. `GFOIL_NOSCATTER` forces the old build (same binary) as
an A/B control.

**Pattern churn (`bench/c1_pattern.py`).** changes/total: golden 1/9, slow cases
8/59…18/56 (max 32% on AH 79-100 C); ≥68% of iterations use the scatter path.

**Gates.** Regression 22/22 bit-identical (scatter and NOSCATTER). 840-case cold
sweep vs post-B baseline: 0 conv / 0 failure_mode / 0 newton_iteration diffs
(identical trajectories). Wall-time (`bench/c4_timing.py`, scatter vs NOSCATTER,
full `run_forward`): golden −3.9%; slow cases −2.2% to −12.8% (GOE 458 −12.8%,
GOE 328 −10.0%), smallest on the highest-churn case (AH 79-100 C −3.2%).

---

## Tape reset A/B: `reset()` replaces `resetHard()` after each AD pass (June 2026)

**What.** Profiling showed 13.9% of AD time in
`std::vector<Direction<3>>::_M_default_append` — the pass-2 adjoint vector being
re-grown on every `run_AD` call because the end-of-pass `tape.resetHard()` in
`ADfuncs.hpp` (`partialOutputspartialInputs`, `partialRpartialx`) released all
chunk/adjoint memory. Both sites now call `tape.reset()`, which clears tape data
and zeroes adjoints but keeps the allocations. The pre-call resets in
`run_AD_py` (`gfoil_ad_bindings.cpp`) are likewise `reset()` — note they were
never the operative ones; A/B-ing them alone changed nothing because the
end-of-pass `resetHard()` had already freed the memory.

**Measurements** (golden case, 30 reps in one process, `bench/tape_reset_ab.py`;
RSS plateau = median of last 10 `VmRSS` readings):

| Variant | AD median (ms) | RSS plateau, 30× run_AD | RSS plateau, 30× fwd+AD alternating |
|---|---|---|---|
| A `resetHard()` | 113.1 | 237.6 MB | — |
| B `reset()` | 86.3 (−24%) | 327.0 MB | 350.1 MB |

**Correctness gate.** Gradients from all 30 B calls are bit-identical to each
other and to all 30 A calls (`reset()` does not leak identifiers/adjoints across
calls); full regression suite 22/22 incl. the forced-transition case.

**Trade.** Each process holds ~90 MB extra steady-state RSS; under a
ProcessPoolExecutor this multiplies by the worker count (e.g. 8 workers ≈
0.7 GB extra), still well under the ~400 MB/process budget. Raw data:
`bench/results/tape_reset_{A,B,B_alt}.json`.

---

## `compute_Dw` as a CoDiPack external function (June 2026)

**What.** `compute_Dw` (`calc_ue_m.hpp`), the wake-influence product
`Dw = Cgam·Bp + Csig` with `Cgam` (Nwake×Ncoords), `Bp` (Ncoords×nPanels),
`Csig` (Nwake×nPanels), `nPanels = Ncoords+Nwake-2`, was being recorded
element-wise on the tape: 30·228·200 = **1,368,000 FMA statements — 31.2% of the
pass-2 tape** (perf profile, June 2026) for a plain matrix product. It is now a
CoDiPack external function following the `implicit_block_solve_b` pattern in the
same file (`ComputeDwData` + `compute_Dw_b` + `compute_Dw_delete`).

**Reverse callback.** Given the output adjoint `L = d̄Dw` (Nwake×nPanels),
looped over `adj->getVectorSize()` dims:

```
d̄Cgam += L · Bpᵀ        (Nwake×Ncoords, Eigen GEMM)
d̄Bp   += Cgamᵀ · L      (Ncoords×nPanels, Eigen GEMM)
d̄Csig += L              (identity, accumulated directly)
```

The EF data stores passive `double` copies of Cgam and Bp (needed in reverse)
plus identifier arrays for Cgam/Bp/Csig (inputs, `updateAdjoint`) and Dw
(outputs, `getAdjoint`+`resetAdjoint`) — ~0.7 MB replacing ~80 MB of tape.

**Primal bit-identity.** The primal is computed with the SAME triple loop in the
SAME summation order as before (i→j, k ascending, then `+ Csig`), on
`.getValue()` doubles — NOT an Eigen GEMM, whose different summation order would
shift Dw at the ~1e-16 level and break the forward golden. Verified: forward
scalars **byte-identical** (out.json `cmp`-equal pre/post). When the tape is
inactive (forward build) the double loop is the entire function, which also
removed compute_Dw's passive-CoDi cost from the forward solve (was 0.6%).
`calc_ue_m`'s row-0 overwrite of Dw (with Bp's last row) is untouched outside
the EF.

**Measured gradient deltas vs the pre-change goldens** (free-transition case):
`dCL/dy` max rel **2.048e-7**, `dCD/dy` **5.168e-8**, `dOASPL/dy` **1.805e-7**;
alpha scalars ~1e-15. Cause: the reverse accumulation for this product now runs
in Eigen GEMM summation order instead of tape statement order — pure
associativity, same class as the 3–6e-8 PanelGeom-refactor precedent. AD-vs-FD
after the change: `dOASPL_dalpha` relerr 5.0e-7 (rtol 1e-11, h=1e-4°) for both
scalar-0.98 and window TE sampling, `dOASPL_dy` at the ~1e-5 FD floor on
O(100)-magnitude nodes (small-|gradient| near-TE nodes show larger *relative*
FD scatter that an h-study confirms is FD noise, not an AD error).

**Wins** (NACA 0012, α=2°, Re=2e6, nCrit=5; same-day before/after, lean -O3
build): pass-2 tape statements 4,382,326 → 3,014,327 (**−31.2%**, exactly the
1.368M compute_Dw statements), Jacobian entries 13.46M → 9.35M, tape memory
258.5 → 178.9 MB (**−79.7 MB**), AD pass wall time 136.8 → 107.8 ms median over
20 reps (**−21.2%**; −19.3% on the 30-iteration mean protocol). Forward wall
time unchanged within noise.

**Goldens regenerated** under this entry's justification: the gradient-array
shifts above are associativity-level and FD-verified, and the forward scalars
are bit-identical. Regenerated together with the regression-suite restoration
below.

---

## Regression suite restored and extended (June 2026)

`tests/` (deleted in bc28bf2, "runtime-only package") is restored from
`bc28bf2~1` and brought up to date: `regression_test.py` now drives
`GFoil.gfoil_cpp` for **two golden cases** — the original free-transition case
and a new **forced-transition case** (same conditions, `transition=[0.1,0.1]`,
golden files `forced_*`) covering the taped-`xift` adjoint path for the first
time (sanity-checked before freezing: `dOASPL_dalpha` AD-vs-FD relerr 3.2e-8 at
rtol 1e-11, h=1e-4°; `dCD_dy` h-study converges to AD, best 1.8e-5 at h=3e-7) —
plus the two **windowed-Amiet OASPL anchors** (scalar 0.98 → 63.18598 dB,
window [0.95,0.99] → 63.36702 dB, tol 1e-6 rel) with their full config embedded
in the test.

**Anchor config clarification.** The "Reference values" paragraph (windowed
Amiet entry, below) omitted the foil: the recorded anchors require the
**analytic open-TE NACA 0012** (4-digit thickness polynomial with the −0.1015
trailing coefficient, cosine spacing, 301 points) — reproduced to ~1e-8.
The `tests/input.json` sharp-TE NACA 0012 coordinates give 63.05499/63.24251
instead (verified identical on a clean 261736a worktree build, so this is foil
geometry, not code drift). The analytic generator is embedded in
`regression_test.py`.

22 checks total; suite passes 22/22 on the current tree, with both golden cases
regenerated here at rel_err exactly 0 by construction (see the compute_Dw entry
above for why regeneration was due).

---

## CMake cleanup + new `noise_run` acoustics-only entry point (June 2026)

Two independent changes.

### 1. Removed the standalone executable builds

`CMakeLists.txt` no longer builds the `GFoil_fwd_codi` (from `src/main.cpp`) and
`GFoil_AD` (from `srcAD/main.cpp`) executables. The Python package
(`pip install .`, scikit-build-core) only ever loads the `gfoil_cpp` pybind11
module and calls `run_forward` / `run_AD`; the two executables were redundant
compiles never used by any Python workflow. The `add_executable(...)` blocks and
all their `target_*` lines were deleted.

`src/main.cpp` and `srcAD/main.cpp` are **not** deleted from disk — they are now
**build-orphaned** (uncompiled) and are candidates for a later `_quarantine/`
sweep per the repo cleanup discipline (not quarantined in this change).

The `gfoil_cpp` module target is unchanged, including its `srcAD/include` and
`srcAD/noise_includes` include paths (still needed by `gfoil_ad_bindings.cpp`),
the `SOURCES` glob, and the three `list(REMOVE_ITEM ...)` lines.

**Verification:** the module configures, compiles, links, and imports; the
forward solve (NACA 0012, α=2°, nCrit=5) and the AD pass (`grad_run`) both run
correctly. No source on the forward/AD path was touched, so the build is
output-identical by construction (the no-aero `noise_run` reproduces the forward
acoustic path bit-for-bit — see below).

### 2. New `noise_run` — acoustics-only, length-generic, never AD'd

A forward-only acoustic investigation entry point that runs the WPS models +
Amiet model directly with **no aerodynamic solve**. The caller supplies
trailing-edge BL states (or a custom wall-pressure spectrum) and gets back the
raw linear wall-pressure spectra and far-field PSD (Pa²/ω — no dB, no
integration, no OASPL). It uses the same `Real = codi::RealReverse` type as the
forward path (so `calc_WPS_*`/Amiet, which need a CoDiPack active type for
`errFunc`, reuse without a `double` retype) but **no gradient is ever taken
through it**.

**Length-generic `_vec` overloads.** The existing fixed-`Nsound` (250) templates
are on the AD-critical forward path and are left **byte-for-byte unchanged**. New
`_vec` overloads accept an arbitrary-length / arbitrary-spacing angular-frequency
array as `std::vector<Real>`:
- `WPSmodels.hpp`: `calc_WPS_{Goody,Kamruzzaman,Rozenburg,Lee,TNO}_vec` —
  identical bodies, `omega[Nsound]`/`phiqq[Nsound]` → `std::vector<Real>`, loop
  bound `Nsound` → `omega.size()`. For TNO the wall-normal arrays stay fixed-size
  (`NblPoints`); only the outer frequency loop becomes runtime-length.
- `sound.hpp`: `calc_WPS_vec` dispatcher — same input floors / `Cf_min` /
  `beta_max` / `H_min` guard logic as `calc_WPS`, dispatching to the `_vec`
  variants.
- `newAmiet.hpp`: `Radiation_integral_total_vec` and `TE_noise_outer_vec` —
  fixed `[Nsound]` arrays/scratch (`C`, `K_bar`, `mu_bar`, `K_1_bar`, `I_abs2`,
  `l_y`) → `std::vector<Real>`; the scalar per-frequency calls
  (`Radiation_integral1/2`, `Estar`, `errFunc`) are reused unchanged.

**Implementation.** `src/include/noise_run.hpp` (`noise_run_cpp<Real>`, template,
forward-only). Inputs: `alphaDeg` (drives the global→TE-local observer rotation,
`te_offset = 0.75·chord`, same transform as `calc_OASPL`), `Re/rho/nu/Ma/chord`
(`Uinf = Re·nu/chord`; the Amiet Mach follows `calc_OASPL`'s `Uinf/340`
convention — `Ma` is accepted but not used for that, so the `_vec` path
reproduces the fixed-size path), observers (global frame, ¼-chord origin), span,
the 7 BL quantities as `[upper, lower]` pairs (`theta, delta_star, tau_max, Ue,
dpdx, tau_wall, delta99`), `freqs_Hz` (any length/spacing), `model`, and an
optional `custom_WPS` (N,2) override.

- **Custom-WPS override:** if given, the BL→WPS path is skipped for both surfaces
  and the columns are used directly (model string irrelevant).
- **Per-surface skip:** BL path skips a surface with `tau_max <= 0` (matches the
  `sound.hpp` gate); custom path skips an all-zero column. A skipped surface has
  a zeros WPS column (no far-field contribution) and passes `Uinf` as its
  fallback edge velocity (mirrors `calc_OASPL`'s `edgeVel = Uinf` fallback).

**Tape hygiene:** `noise_run_cpp` calls `Real::getTape().reset()` at the start to
prevent unbounded tape growth across repeated in-process calls (the forward tape
is inactive by default, so `errFunc`'s `StatementPushHelper` records no
statements — confirmed: 500 repeated calls show zero maxrss growth).

**Bindings / Python:** `run_noise_py` in `gfoil_fwd_bindings.cpp`
(`m.def("noise_run", ...)`); `GFoil.noise_run(...)` wrapper +
`NoiseResult` dataclass (`inputs.py`); both exported from `__init__.py`.

> Note: the `_vec` overloads and `noise_run.hpp` live in `src/include/`
> alongside the existing `WPSmodels.hpp`/`sound.hpp`/`newAmiet.hpp` (the headers
> are there, not in `src/noise_includes/`, which holds only `Faddeeva.hh`).

**Verification.**
- A converged `fwd_run(verbose=True)` (NACA 0012, α=2°, nCrit=5, kam, 2
  observers) fed back into `noise_run` on the **same 250-pt log grid** reproduces
  the forward path's internal `WPS_upper`/`WPS_lower` **and** per-observer
  `FF_spectra` and `obsXYZ_TElocal` to **0.0** difference (exact) — the `_vec`
  path is bit-identical to the fixed-size path on a matching grid.
- Arbitrary length/spacing (37-pt linear grid) runs; custom-WPS override drives
  the far field with BL inputs ignored; an all-zero column zeroes that surface's
  contribution (verified strictly-additive vs the both-surface case at matched
  Ue); a both-zero custom WPS yields zero far-field; `custom_WPS` shape (N,2) is
  validated in the Python wrapper.

---

## Trailing-edge (`x/c == 1.0`) sampling removed (June 2026)

With explicit `[x_lo, x_hi]` windows now available (see *Windowed Amiet TE sample*
below), the legacy `x_target == 1.0` special case in `interpolate_BL_single` — a
hardcoded average over six stations across 0.96–0.985 — is no longer needed and was
removed. `interpolate_BL_single` now always performs the normal single-point
interpolation, which cannot sample the trailing-edge node itself: at `x/c == 1.0`
`find_interp_position` has no node bracket to return (the TE is the boundary), which
is exactly why the averaging hack existed. So sampling at the TE is now **rejected**
rather than silently special-cased:

- Python (`inputs.py`): a scalar `TESampleLoc` must satisfy `0 <= x < 1` (was
  `0 <= x <= 1`); the message points to using a window for TE-region sampling.
- C++ (`extract_BL_TE.hpp` dispatcher): the existing window guard was broadened to
  reject `x_lo >= 1.0 || x_hi >= 1.0` *before* the scalar branch, so both a scalar
  1.0 and a window touching 1.0 throw `std::invalid_argument` (→ Python
  `ValueError`). This covers the standalone/JSON path that bypasses Python.

This **supersedes the earlier Bug-2 `chordScale` fix** (below): that fix corrected
the 1.0 branch, which now no longer exists. Bug 1 (the `get_nodes` shadowing fix)
stands — it is on the normal single-point path.

**Verification** (NACA 0012, α=2°, nCrit=5; baseline = the same pre-window `HEAD`
worktree): scalar `0.98` forward (CL/CD/CM/OASPL) and all four gradient arrays
remain **bit-identical** to `HEAD` (the 0.98 path always took the normal branch, so
deleting the 1.0 branch is a no-op for it); scalar `1.0` cleanly rejected at both
the Python and C++ layers; window `[0.95, 0.99]` unchanged (OASPL bit-matches the
pre-removal windowed build at 63.36701855 @ rtol 1e-11; `dOASPL_dalpha` AD-vs-FD
relerr 5.3e-7, `dOASPL_dy` ~5e-5). No completed optimisation run is affected (all
use 0.98).

---

## Two latent `extract_BL_TE.hpp` bugs fixed (June 2026)

Surfaced during the windowed-Amiet review. Both were **dormant on every completed
run** (all use `TEsample = 0.98`; windows use `x_hi = 0.99`), so the windowed-Amiet
byte-identical scalar verification and gradient gate remain valid — but both were
live landmines for other operating points / chords.

**Bug 1 — variable shadowing in `get_nodes` (bottom surface).** The bottom branch
declared `int botNnodes = botStart;`, a *local* that shadowed the `int& botNnodes`
output reference, so the computed node count was written to a throwaway and the
out-param kept its earlier value (4). Fixed by removing the `int` (assign the
reference, matching the top-surface branch). *Dormant on the 0.98/TE paths because
the bottom stencil there has ≥4 available nodes, so both the buggy value (4) and the
corrected value (clamped to 4) coincide* — confirmed by the scalar-0.98 output and
all four gradient arrays remaining **bit-identical** to the pre-fix `HEAD` build
(`np.array_equal` True). Would have mattered only if a bottom TE stencil offered
<4 turbulent nodes (e.g. transition very near the TE), where the buggy 4 could also
index out of bounds.

**Bug 2 — `chordScale` missing in the `x_target == 1.0` branch.** The normal
single-point path scales `theta`/`delta*` by `chordScale` (they are lengths) before
forming `delta99` and feeding Amiet; the TE special-case averaging branch (the
`NSAMPLES` block) omitted it, so a sample at exactly `x/c = 1.0` fed length scales a
factor `chordScale` wrong into the noise model. Fixed by applying the same
`*= chordScale` to `theta`/`delta*` (before `delta99`) on both surfaces, mirroring
the normal path. **Magnitude:** NACA 0012, α=2°, chord=0.3, `TEsample=1.0`: OASPL
**87.07 dB (buggy) → 81.24 dB (fixed)**, a 5.83 dB error (CL unchanged — only the
acoustic length scales were wrong). Real, but never hit: runs use 0.98, and the
default chord is 1.0 (where `chordScale=1` makes it a no-op anyway). *Superseded:*
the `x_target == 1.0` branch was subsequently removed entirely (see the
TE-sampling-removal entry above), so this fix no longer applies to live code.

**Window/TE interaction guard.** With the 1.0 branch now `chordScale`-consistent, a
*window* must still never let a station land on `x/c == 1.0` — that would route one
station through the nested 0.96–0.985 TE sub-average while the others are point
samples, silently mixing two sampling semantics. The Python layer already rejects
this (`inputs.py` requires `0 < x_lo < x_hi < 1`, so the top station `= x_hi < 1`);
a defensive C++ guard in the window dispatcher now also throws
`std::invalid_argument` (→ Python `ValueError`) if `x_hi >= 1.0`, covering the
standalone/JSON path that bypasses the Python validation. Verified: `[0.96, 1.0]`
rejected at both layers; `[0.95, 0.99]` unchanged (OASPL bit-matches the pre-fix
windowed build, `dOASPL_dalpha` AD-vs-FD relerr 5.3e-7, `dOASPL_dy` ~1e-5).

---

## Windowed Amiet TE sample — BL-averaged wall-pressure input (June 2026)

**Motivation.** A confirmed Kambe–Amiet model exploit: the optimiser drives a
smooth ~0.0022 surface undulation peaked at the *single* wall-pressure sampling
station (`TEsample = 0.98`) to game the predicted OASPL. The undulation is smooth
enough to pass a geometric/curvature constraint (that constraint was calibrated
and correctly rejected — there is no geometric wiggle to catch). The root-cause
fix is model-level: average the Amiet wall-pressure input over a short trailing-
edge window so a localised undulation can no longer swing the predicted OASPL.

**Change.** `Acoustics.TESampleLoc` now accepts *either* a scalar `x/c` (legacy,
default `0.98`) *or* a length-2 `[x_lo, x_hi]` window. For a window the fully
post-processed 7-slot BL/WPS input vectors `[theta, delta*, tau_max, Ue, dpdx,
tau_wall, delta99]` are evaluated at `NWINDOW_SAMPLES` (=9) uniformly spaced
stations across `[x_lo, x_hi]` and **trapezoidally averaged in x/c** (endpoints
interpolated, not snapped), then a single Amiet evaluation runs on the averaged
input — **Option A** (averaged input, one kernel eval), not a distributed-source
integral. The new path reuses the single-point routine verbatim per station, so
each station carries the correct post-processing. (The pre-existing hardcoded
`x_target == 1.0` averaging branch was *not* a clean precedent — it had a latent
`chordScale` bug, fixed separately below.)

**Implementation.** `extract_BL_TE.hpp`: the legacy single-point body was renamed
`interpolate_BL_single()` **unchanged**, and a thin dispatcher
`interpolate_at_95_both_surfaces(..., x_lo, x_hi, ...)` was added. `x_hi <= x_lo`
routes to the single-point routine (byte-identical scalar path); `x_hi > x_lo`
runs the windowed average by calling the same single-point routine per station, so
the two paths cannot drift. A second value `sampleTE_hi` was threaded through
`runCode` (`run_forward.{h,cpp}`), `partialOutputspartialInputs` (`ADfuncs.hpp`),
both pybind layers (`gfoil_fwd_bindings.cpp`, `gfoil_ad_bindings.cpp`, defaulting
to `sampleTE` when the key is absent), and `srcAD/main.cpp`. Python: `inputs.py`
validates a window (`0 < x_lo < x_hi < 1`), `gfoil.py` packs `sampleTE`/
`sampleTE_hi`. The default stays the scalar `0.98`; no opt-script behaviour changed.

**AD/tape hygiene.** `x_lo`/`x_hi` are passive (cast from a double input, never
registered), so the station positions and trapezoid weights are passive constants;
the sampled BL quantities remain taped through `glob.U` and the (geometry-
dependent) node x-positions, so the average is correctly differentiated wrt y and
alpha. All-Real arithmetic — no `fabs`/`hypot`/`.getValue()` introduced.

**Verification** (NACA 0012, α=2°, nCrit=5, Re=2e6, model `kam`; baseline built
from a clean `HEAD` worktree for comparison since the golden/test scaffolding was
dropped in `bc28bf2`):

- *Scalar path byte-identical.* `TEsample=0.98` reproduces baseline `HEAD`
  bit-for-bit: CL/CD/CM/OASPL deltas exactly `0.0`, and the AD gradients
  (`dOASPL_dy`, `dOASPL_dalpha`, `dCL_dy`, `dCL_dalpha`) are bit-identical
  (`np.array_equal` True). The refactor did not perturb the scalar path.
- *Window gradients correct.* AD-vs-central-FD at tightened rtol (1e-11),
  step-size-studied: `dOASPL_dalpha` window relerr **5.3e-7** at `h=1e-4°`
  (scalar control 3.2e-7); `dOASPL_dy` over the TE-window nodes relerr **~1e-5**
  (FD floor, same as the scalar control). Coarser alpha FD steps (≥1e-3°) give
  spurious disagreement from OASPL–alpha curvature — an FD artifact, not a
  gradient bug (confirmed by the h-study converging to AD).
- *Exploit collapse.* On the baseline foil, a smooth Gaussian undulation
  (amp 0.0022, σ=0.004) peaked at x/c=0.98 swings scalar-0.98 OASPL by **+1.06 dB**
  but windowed-[0.95,0.99] OASPL by only **+0.13 dB** (8.4× collapse); an isolated
  single-node bump collapses 13× (+0.294 → +0.022 dB). The window removes the
  single-station lever as intended. (The original ~6 dB figure was the fully-
  optimised w16 runaway, whose artefacts were dropped with the scaffolding; the
  collapse *ratio* is the transferable quantity.) Averaged inputs are bracketed by
  the per-station single-point values by construction (positive trapezoid weights
  summing to 1).

**Reference values** (default rtol 1e-6, observer (1,0,1), span 3): scalar 0.98
OASPL = 63.18598 dB; window [0.95,0.99] OASPL = 63.36702 dB. These are the
window-path regression anchors to add when the test/golden suite is restored
(the scalar anchor is unchanged from `HEAD`).

---

## Forced-transition adjoint fix — `xift` taped (June 2026)

**Bug.** For *forced* transition (`OperatingConds.transition != [1,1]`), the
reverse-mode gradients of CL, CD and OASPL wrt the SVD shape modes disagreed
badly with central finite differences (CD worst: up to ~975% relative error;
CL up to ~9%; OASPL up to ~9%). Free transition was fine.

**Root cause (H1).** The forced-transition arc-length station `xift` is a
continuous function of geometry (via `foil.x` and `distFromStag`), and the BL
residual depends on it (`residuals_shared.hpp`: `xt = param.xift`, weighting the
transition-station interpolation). But `Vsol_t::xift`/`Param_t::xift` were stored
as `double` and `ADfuncs.hpp::partialRpartialx` computed the whole interpolation
with `.getValue()`, so `d(xift)/dy ≡ 0` on the tape. The adjoint term
`−λᵀ ∂R/∂x` was missing `∂R/∂xift · ∂xift/∂x` entirely (absent from both AD
passes). The error was a **flat, step-size-independent floor** — the signature of
a missing constant term — and was largest at the higher (lower-energy) SVD modes
because the missing *absolute* term divided by a smaller gradient. Full diagnosis,
experiments and before/after plots: `bench/results/GRAD_VERIFY.md`.

**Fix.** `Vsol_t::xift` and `Param_t::xift` are now taped `Real`
(`data_structs_shared.hpp`; forward `data_structs.h` stays `double` — value-only).
In `partialRpartialx` the `xift` interpolation is computed in `Real` arithmetic,
moved to **after** `stagpoint_move_AD` so it is measured against
`isol_final.distFromStag` — the *same* arc-length array the residual stations use.
This matches the forward order (`coupled.cpp`:
`stagpoint_move → update_transition → build_glob_RV`) and is essential: the
residual depends on the transition position *relative* to its stations, so an
α-driven stagnation-point shift must cancel between `xift` and the stations. An
intermediate version that used `isol_pre.distFromStag` (inviscid stag) fixed the
geometry modes but **broke `dOASPL/dα`** (0.000% → 3.4%) by leaving a spurious
`dxift/dα`; using `isol_final` restores α to 0.000%. The integer panel-bracket
selection stays passive (`.getValue()` comparisons — legitimately
non-differentiable); only the continuous in-bracket interpolation is taped.

**H2/H3 ruled out for the verification case.** After the H1 fix the OASPL gradient
matches FD to 0.000% for every mode, so no `calc_WPS`/acoustic-floor/`Radiation`
clip is binding and kinking the spectrum here (H2). The α row is 0.000% before and
after, so the radiation-geometry α dependence in `calc_OASPL` (pass 1) is complete
(H3).

**Regression.** `tests/input.json` uses free transition, so the golden is
**bit-identical** before and after this change (10/10, rel_err 0.000e+00); the fix
only affects forced-transition cases. The CLAUDE.md guardrail about a `Real` field
in `Param_t<Real>` did not materialise — `Real xift = 0.0` is a passive constant
until assigned from a taped quantity in forced cases.

---

## API & verbose-output changes (June 2026)

**`ncrithyst` transition-hysteresis logic removed entirely.**
- Decision evidence: the convergence study in `bench/ncrithyst_study/` (seed=0,
  80 aerofoils, 90,720 bare-Newton solves) showed `ncrithyst` does nothing
  meaningful. Convergence is flat within noise at any physically-safe value
  (78.05% at 0 → 78.38% at 0.2). The apparent gains at large values
  (up to 90.5% at 3.0) are an **over-damping artifact** — they pin the
  transition front forward and systematically change the converged solution
  (drag rises monotonically, no plateau), i.e. they buy a convergence number by
  altering the answer, not by stabilising the physics. It was an anti-chatter
  no-op, never a convergence mechanism. See `bench/ncrithyst_study/findings.md`.
- The single hysteresis branch in `update_transition.cpp` (the free-transition
  `else if` arm that pinned `ilam = ilam0` on a spurious single-node retreat,
  with its `amp_first_turb` local) is deleted. Single-node retreats now proceed,
  which is behaviourally **identical to the study's `ncrithyst = 0` ("off")
  column** — the well-characterised intended post-removal behaviour.
- Removed all plumbing: `Param::ncrithyst` (data_structs.h), the `Real ncrithyst`
  parameter + `param.ncrithyst` assignment in `runCode` (run_forward.cpp), the
  `ncrithyst` default in the `runCode` declaration (run_forward.h), the dict read
  and call argument in `gfoil_fwd_bindings.cpp`, the `OperatingConds.ncrithyst`
  field (inputs.py) + `_build_input_dict` entry (gfoil.py), and the dead
  `"ncrithyst"` keys in `tests/input.json` / `input_test.json`. The transient
  `GFOIL_NCRITHYST` constant introduced during the earlier API-removal step is
  also gone.
- AD/adjoint path: `ncrithyst` never appeared in the reverse-mode build
  (`ADfuncs.hpp`, `gfoil_ad_bindings.cpp`); `update_transition.cpp` is a
  forward-only TU. Forward and adjoint transition logic stay consistent and
  `dOASPL/dy`, `dOASPL/dalpha` are unchanged.
- Regression: deleting the branch equals `ncrithyst = 0`, the value the golden
  was generated with, so forward + AD baselines are expected bit-identical.

**Verbose output: per-observer OASPL and TE-local observer coordinates.**
- Two new `VerboseResult` / `ForwardResult` fields, populated only when
  `verbose=True`:
  - `OASPL_perObs` — shape `(nObs,)`, per-observer OASPL [dB re 20µPa].
  - `obsXYZ_TElocal` — shape `(nObs, 3)`, observer coords in the TE-local
    chord-aligned Amiet frame (origin at the trailing edge) [m].
- `OASPL_perObs` integrates the raw linear-power far-field PSD `ff[]`
  (Pa²/(rad/s)) over the linear-spaced ω widths exactly as `calc_OASPL`
  (sound.hpp), including the optional A-weighting and the `1e-30` floor branch.
  Verified: the power-average of `OASPL_perObs` reproduces the published scalar
  `OASPL` to machine precision for both single- and multi-observer cases.

**Acoustic frequency spacing (no change).**
- The acoustic frequency array is, and remains, **log-spaced** (log10-linear
  between `f_min` and `f_max`); see the verbose rebuild in run_forward.cpp and
  the log-spaced grid in calc_OASPL. No code change — recorded here for clarity.

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

### NaN-lock early exit (May 2026) — COMPLETE
Failure mode: at α=±2.5°, nCrit=5, the BL Jacobian goes singular when
the eN ODE reaches near-critical state (amp≈4.91, ilam_top=44) at iter 3.
sparselinsolve NaN fix correctly sets dU=0, but every subsequent iteration
produces the same NaN-triggering Jacobian — frozen BL state means
stagpoint_move and update_transition always produce the same output.
Previously ran all 56 remaining iterations doing nothing before returning
no_convergence.

Fix: added nan_skip_count counter (plain int, no CoDi risk) and
prev_resid_val (plain double) to solve_coupled loop. Detects 3 consecutive
stagnant iterations: residualNorm is NaN, or omega==1.0 with residual not
decreasing by >0.01%. Returns early with failure_mode="nan_lock". Final
failure-mode classifier guarded with !early_exit to prevent overwrite.

Effect: α=±2.5° exits after 7 iterations instead of 59. Total convergence
rate unchanged (cases still need warm-start continuation). NaN entry
confirmed at same fixed Jacobian location every time, confirming frozen state.

Files modified: `src/coupled.cpp`

### init_BL.cpp vsol.forcet[1] clobber fix (June 2026) — COMPLETE
Bug: the wake surface iteration (surf=2) in `init_boundary_layer` wrote
`local_forcet=false` and `local_xift=0.0` into `vsol.forcet[1]` and
`vsol.xift[1]`, overwriting the upper-surface forced-transition state set
in the surf=1 pass. The `if (surf < 2 && ...)` guard meant `local_forcet`
and `local_xift` were never updated for the wake, so surf=2 always clobbered
surf=1's values with false/0.

Effect: `build_glob_RV` (called before `update_transition` on iteration 1)
read the clobbered `vsol.forcet[1]=false` and assembled the first Jacobian
with the free-transition path instead of forced. `update_transition` repaired
the value from iteration 2 onwards, so convergence was not broken but the
first Newton Jacobian was wrong when forced transition was active.

Fix: guard the write with `if (surf < 2)` so the wake iteration never touches
`vsol.forcet` or `vsol.xift`.

Files modified: `src/init_BL.cpp`

---

## Known Limitations (full write-ups)

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

### Additional multi-node BL attractor cases — NACA 0012, nCrit=5 (May 2026) — ACCEPTED LIMITATIONS

Two further cold-start failure modes identified and characterised on NACA 0012,
nCrit=5 during investigation of newly-backstepping angles.

**α=±2.5° — NaN-lock (fast exit implemented):**
The Newton path reaches ilam_top=44 with amp≈4.91 at iteration 3, placing the
eN ODE in a near-critical state that produces a singular BL Jacobian block
(confirmed at fixed entries: rows 180/417, cols 244/552). The sparselinsolve
NaN fix sets dU=0, but stagpoint_move and update_transition see the same frozen
BL state each time and produce the same output — every subsequent iteration
re-triggers the same singular Jacobian. Previously burned all 56 remaining
iterations returning no_convergence. Now returns failure_mode="nan_lock" after
7 iterations (see NaN-lock early exit fix). Warm-start from a neighbouring
alpha converges without issue. Most likely cause: the period-2 fix in
update_transition.cpp (ilam==ilam0 branch) shifted the iter-3 BL state to
this NaN-triggering configuration.

**α=4.7° — multi-node BL attractor (same class as NACA 0008-34):**
Distinct from the single-node period-2 oscillation. omega drops to 0.001 by
iteration 7 with ilam stable from iteration 1 — the ctau limiter is saturated
before the freeze fires. Both surfaces independently satisfy the primary freeze
condition at iteration 10 (stable_ilam_iters=9 ≥ 8, and current residual 0.171
exceeds iter-2 oldest/2 = 0.141). Co-activation plays no role — both surfaces
freeze simultaneously via the primary condition regardless of the co-activation
threshold. Confirmed by disabling co-activation entirely: identical failure.
The freeze damps ctau oscillations but the residual plateaus at ~0.17 (340×
larger than NACA 0008-34 which reached ~5×10⁻⁴). Root cause: the Newton step
direction at this BL state produces large ctau corrections that the 0.05/step
limiter caps every iteration, preventing net progress. Warm-start from
α=4.5° or α=5.0° converges correctly. failure_mode="transition_front_oscillation"
is returned; backstepping continuation handles it automatically.

Do NOT attempt to fix α=4.7° via the co-activation threshold — the two surfaces
freeze simultaneously through the primary condition, making the co-activation
threshold irrelevant.

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

## Deduplication History (May 2026)

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

---

## Feature Implementation Notes

### Multiple observer locations (COMPLETE)
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

### A-weighting toggle (COMPLETE)
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

### Warm-start fix in pybind11 path (COMPLETE)
- runCode gains `const RestartState* warmStart = nullptr` parameter
- When non-null, initialises glob.U and vsol.turb from in-memory state
  instead of reading restart.json
- run_forward_py gains optional `prev_jacobian` argument
- _call_forward in gfoil.py passes prev_result states/turb through
- standard_run tracks last_converged and passes it as warm start
  in the forward-stepping loop
- Subprocess fallback: writes restart.json from prev_result if available
- Binary path unchanged: fromRestart=1 still reads restart.json

### Pybind11 Python bindings (COMPLETE)
- `gfoil_cpp.cpython-38-*.so` built into `GFoil/GFoil/` (importable as `from . import gfoil_cpp`)
- `fwd_run()` returns `FwdResult`; `grad_run(result, ...)` returns `GradResult`
- No file I/O in pybind11 path
- Build: `cmake -B build . -DPYBIND11_PYTHON_VERSION=3.8 -DPYTHON_EXECUTABLE=$(pyenv which python3.8) && cmake --build build --target gfoil_cpp`

### newAmiet.hpp physics fixes (supervisor review)
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
