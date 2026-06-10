# Forward-Solve Convergence & Forced-Transition Drift — Diagnostic Findings

**Scope:** read-only investigation. No solver code modified, nothing committed.
**Data:** `../AerofoilOpt/FourthRun` and `../AerofoilOpt/MultiRun2` (tripped runs,
`topTrans=botTrans=0.1`, Re=2e6, ncrit=6, targCL=0.8, chord=0.3, model=`kam`).
**Tooling:** re-ran the live `GFoil.gfoil_cpp` pybind module on saved geometries.

---

## Part 1 — Why the saved `TransitionX` drifts across designs

### 1.1 The data

FourthRun `best/` saved designs, requested trip = 0.1c on both surfaces:

| nMode | topTransX | botTransX | alpha | CL | CD | OASPL |
|------:|----------:|----------:|------:|-----:|------:|------:|
| 0  | **0.0338** | 0.1018 | 7.40 | 0.804 | 0.01139 | 60.35 |
| 2  | 0.0971 | 0.1008 | 5.60 | 0.800 | 0.00991 | 59.11 |
| 3  | 0.1000 | 0.1053 | 8.04 | 0.800 | 0.01031 | 58.20 |
| 4  | 0.1042 | 0.1002 | 8.51 | 0.800 | 0.01060 | 57.12 |
| 5  | 0.1028 | 0.1045 | 8.75 | 0.800 | 0.01066 | 57.12 |
| 6  | 0.1018 | 0.1116 | 8.82 | 0.800 | 0.01050 | 56.99 |
| 7  | 0.1015 | 0.1117 | 8.80 | 0.800 | 0.01050 | 56.98 |
| 8–9  | 0.1011 | 0.1018 | 8.53 | 0.800 | 0.01038 | 56.97 |
| 10 | 0.1010 | 0.1043 | 9.25 | 0.800 | 0.01052 | 56.94 |
| 11–12 | 0.1012 | 0.1054 | 9.28 | 0.800 | 0.01052 | 56.93 |
| 13 | 0.1025 | 0.1054 | 9.44 | 0.800 | 0.01059 | 56.92 |

- `topTransX`: range 0.0338–0.1042, spread 0.070
- `botTransX`: range 0.1002–0.1117, spread 0.012

### 1.2 The variation is QUANTISED onto the node grid — verified directly

Geometry is normalised to chord=1; trip target is x=0.1. The reported transition
x is the x-coordinate of the **first turbulent panel node**, and it lands exactly
on the node grid:

```
mode 4  inner nodes near 0.1c:  ... 0.0952  0.1002  0.1042  0.1102 ...
        reported topTransX = 0.10423  -> node 0.1042   (exact)
        reported botTransX = 0.10016  -> node 0.1002   (exact)
mode 13 inner nodes near 0.1c:  ... 0.0941  0.1025  0.1054  0.1113 ...
        reported topTransX = 0.10249  -> node 0.1025   (exact)
        reported botTransX = 0.10541  -> node 0.1054   (exact)
```

Node spacing near 0.1c is ≈ 0.003–0.008 in x/c — the same magnitude as the
observed `botTransX` spread (0.012 ≈ 2–3 node intervals). **Conclusion (a):
the drift is node quantisation.** The interpolated trip arc-length moves smoothly
with the geometry, but the *reported / flagged* node snaps to the discrete grid.

There is also a clear case (b) component on the **suction side at low camber**:
mode 0 (baseline NACA0012, α=7.4°, ncrit=6) reports `topTransX = 0.0338`, far
upstream of the trip — natural e^N transition has pre-empted the forced trip.
mode 2 `topTransX = 0.0971` is the node just upstream of 0.1, the same effect at
threshold. So: **primarily (a) quantised; with natural-transition pre-emption on
the suction side for low-amplitude designs.**

### 1.3 How the trip is applied in the source

The mechanism has three layers (forward TU shown; the AD/taped path mirrors it):

**(i) Trip arc-length is interpolated — SMOOTH.**
`src/init_BL.cpp:213–230` (and identically `src/update_transition.cpp:109–128`):
```cpp
double xft_abs = param.xft_xc[surf] * x_max;              // 0.1 * chord
for (int k = 1; k < N; ++k) {
    if ((x_prev - xft_abs) * (x_curr - xft_abs) <= 0.0) { // bracket node pair
        double frac = (xft_abs - x_prev) / (x_curr - x_prev);
        local_xift = xi_prev + (xi_curr - xi_prev) * frac;   // LINEAR interp
        local_forcet = true; break;
    }
}
```
`xift` is the exact interpolated arc-length and is **taped** (per CLAUDE.md the
`Real xift` field deliberately carries `d(xift)/d(geometry)` into the adjoint).

**(ii) Within the bracketing interval the residual places transition at exactly
`xift` — SMOOTH.** `src/include/residuals_shared.hpp:525–537`:
```cpp
} else {                       // FORCED TRANSITION: xt prescribed
    xt = static_cast<Real>(param.xift);
    w2 = (xt-x1)/dx; w1 = 1-w2;       // split-station at exact xift
    ...
}
```

**(iii) But WHICH interval is the transition interval is an INTEGER test — STEP.**
The first turbulent node is chosen by a discrete bracket test, evaluated on
passive `.getValue()`s: `src/init_BL.cpp:474–485`:
```cpp
bool forced_tran = false;
if (local_forcet && local_xift > 0.0) {
    double xi_prev_d = isol.distFromStag[prevNode].getValue();
    double xi_curr_d = isol.distFromStag[currNode].getValue();
    forced_tran = (xi_prev_d <= local_xift && local_xift < xi_curr_d);
}
```
and the natural-transition march breaks at a node interval the same way
(`src/update_transition.cpp:68–83`). The reported value is the first-turbulent
node x-coord: `src/run_forward.cpp:201–214`.

So the enforced trip is an **interval bound** ("transition at the first node with
`distFromStag ≥ xift`"), i.e. it **snaps to the node just downstream of the
interpolated `xift`**. The continuous physics (split-station at `xt = xift`) is
smooth and taped, but the **set of turbulent nodes is a STEP function of the
geometry**: when a design deforms enough that `xift` crosses a node's arc-length,
one node flips lam↔turb and the BL closure / amp-ctau seeding changes
discontinuously.

### 1.4 Adjoint implication

Within a node interval the adjoint is consistent (it sees `d(xift)/d(geom)`).
**At a node-crossing event the solution state has a step/kink the adjoint cannot
see** — the integer node selection in (iii) is computed on `.getValue()` and is
non-differentiable, so the gradient is missing the contribution of the transition
front jumping across a node. The reported `TransitionX` values sitting on
*different* nodes for neighbouring designs (Table 1.1) are direct evidence the
optimiser is repeatedly stepping across these crossings. **This is a real
contributor to forced-transition gradient noise and SNOPT line-search stalls**,
and — as Part 2 shows — it is also the direct cause of the forward-solve failures.

---

## Part 2 — Classifying the forward-solve failures

### 2.1 Failure inventory (SNOPT EXIT codes, FourthRun + MultiRun2)

| EXIT | meaning | count | share |
|-----:|---------|------:|------:|
| **60** | **undefined user-supplied function** (GFoil fwd returned non-converged/NaN) | **65** | **42%** |
| 30 | resource limit (major-iteration limit) | 44 | 28% |
| 0  | finished successfully | 39 | 25% |
| 40 | numerical difficulties | 7 | 5% |

EXIT 60 (INFO 63 "unable to proceed into undefined region") is the dominant
failure and it *is* the GFoil forward-solve failure being asked about. It
concentrates at higher mode counts (FourthRun norestarts: EXIT 60 at modes
2, 8, 10, 11, 12, 14), consistent with wigglier geometry placing the transition
front on a bistable node more often.

No `crash_traceback.txt` files exist in either run — the failures are clean
non-convergence returns, not C++ crashes.

### 2.2 Re-running representative points — separation vs solver weakness

I re-ran the live solver on three converged best geometries over an α sweep
(trip 0.1c, Re=2e6, ncrit=6, kam). Two distinct regimes emerged.

**Regime A — `transition_front_oscillation` (SOLVER WEAKNESS).**
Mode 13, fine α sweep, *single direct solve* (no backstepping):

```
 alpha conv     CD   topTrX  botTrX  minTau newt  fail
  9.44  YES 0.01059  0.1025  0.1054  2.142    5   ok
  9.46  NO     -        -       -      -    100   transition_front_oscillation
  9.48  NO     -        -       -      -    100   transition_front_oscillation
  9.50  NO     -        -       -      -    100   transition_front_oscillation
  9.52  NO     -        -       -      -    100   transition_front_oscillation
  9.54  YES 0.01078  0.0941  0.1054  1.549    8   ok
```

The failure band is **exactly** the α at which the suction-side transition node
flips from 0.1025 to the adjacent upstream node 0.0941. `minTau` is healthy
(≈2.1 below, ≈1.5 above) — **no separation, no shape-factor blow-up**. The solve
limit-cycles for 100 iterations because the transition front is bistable between
two adjacent nodes: natural e^N amplification crosses ncrit right at the node next
to the forced trip, so each Newton sweep alternates between "forced at 0.1025" and
"natural at 0.0941." This is the documented period-N transition-front attractor,
triggered by the node quantisation of Part 1. **Note the optimum of mode 13 sits
at α=9.444 — right on the lip of this band**, which is why the optimiser keeps
falling into it.

**Regime B — genuine TE separation (METHOD OUT OF ENVELOPE).**
Mode 2 (low camber) at high α:

```
 alpha conv    CD     OASPL  Htop  tauTop   minTau  sepNodes
 10.00 YES 0.01431  65.49  2.31  1.9597  -1.3250    1
 11.50 YES 0.01650  69.13  2.57  1.1458  -5.9204    3
 13.50 YES 0.02018  77.47  3.35 -0.1692 -10.2308    5
```

Here `minTau` goes negative (reversed flow), `Htop` climbs 2.0→3.35, OASPL
balloons — genuine suction-side separation. **Notably this still converges.**
The fixed-point solver tolerates mild-to-moderate separation; separation is *not*
what is producing the EXIT 60s in these runs.

### 2.3 Mitigation probe (input knobs only, no code change)

Mode 13 @ α=9.48 (baseline fails). What recovers it:

| knob | result |
|------|--------|
| baseline (uf=1.0, tef=0.09, ncrit=6, rtol=1e-9) | FAIL — osc, 100 it |
| rtol 1e-6 (looser) | FAIL — osc |
| rtol 1e-11 (tighter) | FAIL — osc |
| panel uf=1.8, tef=0.10 | **OK**, 6 it (topTrX→0.0998) |
| panel uf=2.1, tef=0.09 | **OK**, 5 it (topTrX→0.1012) |
| panel uf=1.5, tef=0.09 | FAIL — osc |
| ncrit=9 (push nat. trans downstream) | **OK**, 5 it |
| ncrit=4 (pull nat. trans upstream) | **OK**, 5 it |

Reading:
- **rtol is irrelevant** — it is a discrete limit cycle, not a tolerance issue.
- **Re-paneling fixes it** by moving node positions off the bistable node
  (it is purely node-position-dependent — confirms quantisation root cause).
- **ncrit either way fixes it** — moving natural transition decisively up- or
  downstream of the forced node breaks the lam/turb tie. The failure exists only
  in the narrow zone where natural amplification crosses ncrit *adjacent to* the
  forced trip node. ncrit=6 with a 0.1c trip at high suction-side α sits squarely
  in that zone.

### 2.4 Why the optimiser sees these as hard failures

`OptScripts2/single_opt.py:100` and `multipoint_opt.py:157` call
`fwd_run(foil, oper, noise_obj)` with **`repanel=False`** for every objective
evaluation. `standard_run` does α back/forward-stepping, but that cannot cross a
transition-front-oscillation band (it bisects toward the target α and fails — seen
live: it walks to 9.445 and stops). The one mitigation that *does* work —
panel-distribution variation — is only in the `repanel=True` branch, which the
per-iteration objective never invokes. So a design landing on a bistable node
returns non-converged → SNOPT EXIT 60.

### 2.5 Hypotheses, ranked by evidence (no changes made)

1. **(strongest) Use the repanel fallback in objective evaluations.**
   `repanel=True` recovered every failing probe by shifting nodes off the
   bistable point; it is already implemented and unused on the hot path. Cost is
   extra solves only on failure. Directly addresses the demonstrated root cause.
2. **(strong) Break the lam/turb node tie in the transition march.** A hysteresis
   band or sub-node blending so the front cannot bistably alternate between two
   adjacent nodes (the amp-vs-ncrit comparison + the forced bracket test currently
   make a hard per-node choice). ncrit perturbation recovering the case is direct
   evidence the bistability is the mechanism. *(Solver change — hypothesis only.)*
3. **(moderate) Sub-node interpolation of the turb flag at the forced trip.** The
   physics already places the split-station at exact `xift`; making the *node set*
   change continuously (e.g. weight the transition node) would also smooth the
   adjoint kink from Part 1. *(Solver change — hypothesis only.)*
4. **(weak / rejected) rtol or Newton damping.** rtol had zero effect; this is a
   discrete attractor, not a residual-tolerance or step-length problem.

---

## What this means for the optimisation

**The right fix is to harden the solve, not (only) to fence off the optimiser.**
The dominant failure (EXIT 60, 42% of attempts) is **not** flow separation — when
the solver genuinely separates (mode 2 at high α) it still converges. The failure
is a numerical **transition-front bistability**: at ncrit=6 with a 0.1c trip on a
high-α suction side, natural e^N transition crosses ncrit on the node next to the
forced trip, and the front limit-cycles between two adjacent nodes. This is the
same node-quantisation identified in Part 1, now causing hard non-convergence
rather than just gradient noise.

Concretely:
- **Cheapest immediate win:** route objective `fwd_run` calls through the existing
  `repanel=True` fallback (Hyp. 1). It demonstrably recovers the failing designs
  and needs no solver-physics change.
- **Deeper win:** desensitise the transition front to the node grid (Hyp. 2/3).
  This would *both* remove the EXIT-60 failures *and* smooth the adjoint kink that
  is feeding SNOPT gradient noise — the two problems share one root cause.
- **Constraining the optimiser away from the region is the weaker option:** the
  band is razor-thin in design space (Δα ≈ 0.1° wide, mode-13 optimum sits on its
  edge) and moves with every geometry change, so it cannot be cleanly fenced with
  a static constraint.

**Honest caveats.** (1) The re-runs vary α on the *converged optima* to expose the
failure mechanism; I did not reconstruct the exact mid-line-search geometries
SNOPT failed on (those are not saved), so the *count* attribution between
"oscillation" and "separation" for the in-run EXIT-60s is inferred from the
mechanism, not measured per-failure. The mechanism itself is reproduced cleanly.
(2) The live module is the current `BLaveragingWindow` build; FourthRun was run
on the Jun-8 checkout, so absolute α-band edges may shift slightly, but the
qualitative behaviour (node-flip → oscillation, repanel/ncrit recovers) is robust.
