# Cold-start convergence baseline — 60-foil sweep

Library: `Smoothed_TEfixed_linear/` (1300 files); 60 foils stratified-sampled (seed=0, sorted order, one per contiguous stratum).
Grid: alpha [-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0] x nCrit [5, 9], Re=2e+06, rtol=1e-06, repanel=False, free transition, model kam, TESampleLoc 0.98, observer (0,3,0.5), span 3.
Cold start = single `gfoil_cpp.run_forward` call (what fwd_run tries first); `rescue_conv` = full fwd_run backstepping outcome for cold failures.

Rerun: `python3 bench/cold_start_sweep.py` (venv: /home/pa20830/NewGradientVal/venv)

## Headline

- cold-converged (valid OASPL): **707/840 = 84.2%**
- true aero cold-converged (conv + acoustic_nan): **707/840 = 84.2%**
- cold failures rescued by fwd_run backstepping: 97/133
- Newton iterations (converged): median 16, p90 41, max 59
- wall ms (converged): median 151, p90 351, max 635

## Convergence by alpha

| alpha | cold conv | true aero |
|---|---|---|
| -4 | 83.3% | 83.3% |
| -2 | 85.8% | 85.8% |
| +0 | 91.7% | 91.7% |
| +2 | 92.5% | 92.5% |
| +4 | 88.3% | 88.3% |
| +6 | 78.3% | 78.3% |
| +8 | 69.2% | 69.2% |

## Convergence by nCrit

| nCrit | cold conv | true aero |
|---|---|---|
| 5 | 84.0% | 84.0% |
| 9 | 84.3% | 84.3% |

## Failure modes

| failure_mode | count | % of grid | rescued by backstepping |
|---|---|---|---|
| transition_front_oscillation | 54 | 6.4% | 40/54 |
| nan_lock | 32 | 3.8% | 16/32 |
| no_convergence | 27 | 3.2% | 23/27 |
| diverged | 20 | 2.4% | 18/20 |

## Slowest-converging 10 cases

| foil | alpha | nCrit | iters | ms |
|---|---|---|---|---|
| S4083 (8%).dat | +6 | 9 | 59 | 457 |
| BELL-WORTMANN FX 69-H-083 AIRFOIL.dat | +4 | 5 | 58 | 464 |
| GOE 303 (FRIEDRICHSHAFEN G03) AIRFOIL.dat | +0 | 9 | 58 | 551 |
| GOE 328 AIRFOIL.dat | +8 | 9 | 58 | 482 |
| GOE 458 AIRFOIL.dat | -2 | 9 | 57 | 535 |
| NASA-LANGLEY MS(1)-0313 AIRFOIL.dat | +8 | 9 | 57 | 436 |
| WORTMANN FX 63-137 AIRFOIL.dat | +6 | 9 | 57 | 417 |
| AH 79-100 C AIRFOIL.dat | +6 | 9 | 56 | 448 |
| MH 93  15.98%.dat | -4 | 5 | 56 | 480 |
| E374.dat | -2 | 9 | 55 | 437 |

## All failing cases (file, alpha, nCrit, mode, rescued)

- `AH 79-100 C AIRFOIL.dat` alpha=+8 nCrit=9 — diverged (rescued)
- `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` alpha=+6 nCrit=9 — diverged (rescued)
- `CH10 (smoothed).dat` alpha=+8 nCrit=5 — diverged (rescued)
- `EPPLER 1098 AIRFOIL.dat` alpha=+0 nCrit=5 — diverged (rescued)
- `FX 61-168 AIRFOIL.dat` alpha=+8 nCrit=5 — diverged (rescued)
- `FX60-100 10.0% smoothed.dat` alpha=+8 nCrit=9 — diverged (rescued)
- `GOE 117 (MVA MK.4) AIRFOIL.dat` alpha=+6 nCrit=9 — diverged
- `GOE 184 (MVA H.29) AIRFOIL.dat` alpha=-4 nCrit=9 — diverged (rescued)
- `GOE 346 (FRIEDRICHSHAFEN-STAAKEN) AIRFOIL.dat` alpha=+8 nCrit=5 — diverged (rescued)
- `GOE 393 AIRFOIL.dat` alpha=-4 nCrit=9 — diverged (rescued)
- `GOE 458 AIRFOIL.dat` alpha=+4 nCrit=9 — diverged (rescued)
- `GOE 513 AIRFOIL.dat` alpha=+8 nCrit=9 — diverged (rescued)
- `GOE 596 AIRFOIL.dat` alpha=+8 nCrit=9 — diverged (rescued)
- `GOE 777 AIRFOIL.dat` alpha=+8 nCrit=9 — diverged (rescued)
- `LWK 80-150-K25.dat` alpha=+4 nCrit=5 — diverged
- `MH 121  8.76%.dat` alpha=-2 nCrit=5 — diverged (rescued)
- `MH 121  8.76%.dat` alpha=+2 nCrit=9 — diverged (rescued)
- `NASA SC(2)-0012 AIRFOIL.dat` alpha=-4 nCrit=5 — diverged (rescued)
- `NASA SC(2)-0012 AIRFOIL.dat` alpha=+4 nCrit=5 — diverged (rescued)
- `S9033 (7.5%).dat` alpha=+8 nCrit=5 — diverged (rescued)
- `AH 79-100 C AIRFOIL.dat` alpha=+6 nCrit=5 — nan_lock (rescued)
- `CH10 (smoothed).dat` alpha=+8 nCrit=9 — nan_lock (rescued)
- `E374.dat` alpha=+4 nCrit=9 — nan_lock (rescued)
- `E374.dat` alpha=+8 nCrit=5 — nan_lock (rescued)
- `EPPLER 502 AIRFOIL.dat` alpha=-4 nCrit=9 — nan_lock (rescued)
- `EPPLER 502 AIRFOIL.dat` alpha=+8 nCrit=5 — nan_lock (rescued)
- `EPPLER 864 STRUT AIRFOIL.dat` alpha=-2 nCrit=5 — nan_lock (rescued)
- `EPPLER 864 STRUT AIRFOIL.dat` alpha=+6 nCrit=9 — nan_lock
- `FX 61-168 AIRFOIL.dat` alpha=-4 nCrit=9 — nan_lock (rescued)
- `FX60-100 10.0% smoothed.dat` alpha=-2 nCrit=5 — nan_lock (rescued)
- `GOE 443 AIRFOIL.dat` alpha=-4 nCrit=5 — nan_lock (rescued)
- `GOE 443 AIRFOIL.dat` alpha=+4 nCrit=5 — nan_lock (rescued)
- `GOE 513 AIRFOIL.dat` alpha=+6 nCrit=9 — nan_lock (rescued)
- `GOE 777 AIRFOIL.dat` alpha=+8 nCrit=5 — nan_lock (rescued)
- `MH 121  8.76%.dat` alpha=+2 nCrit=5 — nan_lock
- `NACA 2415.dat` alpha=+8 nCrit=5 — nan_lock (rescued)
- `OAF128 AIRFOIL.dat` alpha=-4 nCrit=5 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=-4 nCrit=9 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=-2 nCrit=5 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=-2 nCrit=9 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+0 nCrit=5 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+0 nCrit=9 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+2 nCrit=5 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+2 nCrit=9 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+4 nCrit=5 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+4 nCrit=9 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+6 nCrit=5 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+6 nCrit=9 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+8 nCrit=5 — nan_lock
- `OAF128 AIRFOIL.dat` alpha=+8 nCrit=9 — nan_lock
- `SD6080 (9.2%).dat` alpha=+0 nCrit=5 — nan_lock (rescued)
- `USA 5 AIRFOIL.dat` alpha=+6 nCrit=9 — nan_lock (rescued)
- `AH 79-100 C AIRFOIL.dat` alpha=+0 nCrit=9 — no_convergence (rescued)
- `CH10 (smoothed).dat` alpha=-2 nCrit=9 — no_convergence (rescued)
- `EPPLER 1098 AIRFOIL.dat` alpha=+8 nCrit=5 — no_convergence (rescued)
- `EPPLER 502 AIRFOIL.dat` alpha=+6 nCrit=9 — no_convergence (rescued)
- `EPPLER 582 AIRFOIL.dat` alpha=+8 nCrit=9 — no_convergence (rescued)
- `EPPLER 864 STRUT AIRFOIL.dat` alpha=-2 nCrit=9 — no_convergence
- `EPPLER 864 STRUT AIRFOIL.dat` alpha=+0 nCrit=9 — no_convergence
- `EPPLER 864 STRUT AIRFOIL.dat` alpha=+4 nCrit=5 — no_convergence (rescued)
- `FX 61-168 AIRFOIL.dat` alpha=-4 nCrit=5 — no_convergence (rescued)
- `FX 67-K-150-17 AIRFOIL.dat` alpha=+6 nCrit=5 — no_convergence (rescued)
- `FX 67-K-150-17 AIRFOIL.dat` alpha=+8 nCrit=9 — no_convergence
- `FX60-100 10.0% smoothed.dat` alpha=+6 nCrit=9 — no_convergence (rescued)
- `GOE 226 (MVA H.36) AIRFOIL.dat` alpha=-2 nCrit=9 — no_convergence (rescued)
- `GOE 513 AIRFOIL.dat` alpha=-4 nCrit=5 — no_convergence (rescued)
- `GOE 513 AIRFOIL.dat` alpha=-4 nCrit=9 — no_convergence (rescued)
- `GOE 602 MOD. AIRFOIL.dat` alpha=+6 nCrit=9 — no_convergence (rescued)
- `GOE 629 AIRFOIL.dat` alpha=-2 nCrit=5 — no_convergence (rescued)
- `GOE 629 AIRFOIL.dat` alpha=+8 nCrit=9 — no_convergence (rescued)
- `GOE 770 AIRFOIL.dat` alpha=+6 nCrit=9 — no_convergence (rescued)
- `GOE 777 AIRFOIL.dat` alpha=+0 nCrit=9 — no_convergence (rescued)
- `HQ 2.5-8 AIRFOIL.dat` alpha=+6 nCrit=5 — no_convergence (rescued)
- `LWK 80-150-K25.dat` alpha=-4 nCrit=5 — no_convergence
- `MH 121  8.76%.dat` alpha=-2 nCrit=9 — no_convergence (rescued)
- `NACA 2415.dat` alpha=-4 nCrit=9 — no_convergence (rescued)
- `NACA 65-210.dat` alpha=+8 nCrit=9 — no_convergence (rescued)
- `S4083 (8%).dat` alpha=+8 nCrit=5 — no_convergence (rescued)
- `WORTMANN FX 63-137 AIRFOIL.dat` alpha=+8 nCrit=9 — no_convergence (rescued)
- `AG17.dat` alpha=+6 nCrit=5 — transition_front_oscillation (rescued)
- `AG17.dat` alpha=+6 nCrit=9 — transition_front_oscillation (rescued)
- `AH 79-100 C AIRFOIL.dat` alpha=-2 nCrit=9 — transition_front_oscillation (rescued)
- `AH21 7% version (Andrew Hollom).dat` alpha=-4 nCrit=5 — transition_front_oscillation (rescued)
- `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` alpha=+0 nCrit=5 — transition_front_oscillation (rescued)
- `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` alpha=+2 nCrit=9 — transition_front_oscillation (rescued)
- `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` alpha=+6 nCrit=5 — transition_front_oscillation (rescued)
- `E186 (10.27%).dat` alpha=-2 nCrit=5 — transition_front_oscillation
- `EPPLER 1098 AIRFOIL.dat` alpha=-4 nCrit=5 — transition_front_oscillation (rescued)
- `EPPLER 399 AIRFOIL.dat` alpha=+6 nCrit=5 — transition_front_oscillation
- `EPPLER 399 AIRFOIL.dat` alpha=+8 nCrit=5 — transition_front_oscillation
- `EPPLER 399 AIRFOIL.dat` alpha=+8 nCrit=9 — transition_front_oscillation (rescued)
- `EPPLER 502 AIRFOIL.dat` alpha=+4 nCrit=9 — transition_front_oscillation (rescued)
- `EPPLER 502 AIRFOIL.dat` alpha=+8 nCrit=9 — transition_front_oscillation (rescued)
- `EPPLER 637 AIRFOIL.dat` alpha=+0 nCrit=5 — transition_front_oscillation (rescued)
- `EPPLER 864 STRUT AIRFOIL.dat` alpha=+2 nCrit=9 — transition_front_oscillation
- `EPPLER 864 STRUT AIRFOIL.dat` alpha=+4 nCrit=9 — transition_front_oscillation
- `FX 61-168 AIRFOIL.dat` alpha=+8 nCrit=9 — transition_front_oscillation (rescued)
- `FX 67-K-150-17 AIRFOIL.dat` alpha=+8 nCrit=5 — transition_front_oscillation (rescued)
- `FX60-100 10.0% smoothed.dat` alpha=+4 nCrit=5 — transition_front_oscillation (rescued)
- `FX60-100 10.0% smoothed.dat` alpha=+4 nCrit=9 — transition_front_oscillation (rescued)
- `FX60-100 10.0% smoothed.dat` alpha=+8 nCrit=5 — transition_front_oscillation (rescued)
- `GOE 117 (MVA MK.4) AIRFOIL.dat` alpha=+6 nCrit=5 — transition_front_oscillation (rescued)
- `GOE 117 (MVA MK.4) AIRFOIL.dat` alpha=+8 nCrit=5 — transition_front_oscillation (rescued)
- `GOE 328 AIRFOIL.dat` alpha=+2 nCrit=5 — transition_front_oscillation
- `GOE 346 (FRIEDRICHSHAFEN-STAAKEN) AIRFOIL.dat` alpha=+8 nCrit=9 — transition_front_oscillation (rescued)
- `GOE 393 AIRFOIL.dat` alpha=-4 nCrit=5 — transition_front_oscillation (rescued)
- `GOE 425 AIRFOIL.dat` alpha=+0 nCrit=5 — transition_front_oscillation (rescued)
- `GOE 443 AIRFOIL.dat` alpha=+6 nCrit=5 — transition_front_oscillation
- `GOE 443 AIRFOIL.dat` alpha=+8 nCrit=5 — transition_front_oscillation
- `GOE 770 AIRFOIL.dat` alpha=-4 nCrit=5 — transition_front_oscillation (rescued)
- `GOE 770 AIRFOIL.dat` alpha=-2 nCrit=5 — transition_front_oscillation (rescued)
- `GOE 777 AIRFOIL.dat` alpha=-4 nCrit=9 — transition_front_oscillation (rescued)
- `GOE 777 AIRFOIL.dat` alpha=-2 nCrit=5 — transition_front_oscillation (rescued)
- `HQ 1.5-8 AIRFOIL.dat` alpha=+2 nCrit=5 — transition_front_oscillation (rescued)
- `HQ 1.5-8 AIRFOIL.dat` alpha=+2 nCrit=9 — transition_front_oscillation (rescued)
- `HQ 1.5-8 AIRFOIL.dat` alpha=+4 nCrit=5 — transition_front_oscillation
- `HQ 1.5-8 AIRFOIL.dat` alpha=+6 nCrit=9 — transition_front_oscillation (rescued)
- `HQ 2.5-8 AIRFOIL.dat` alpha=-2 nCrit=9 — transition_front_oscillation (rescued)
- `LWK 80-150-K25.dat` alpha=+8 nCrit=5 — transition_front_oscillation (rescued)
- `MH 49.dat` alpha=+8 nCrit=9 — transition_front_oscillation (rescued)
- `NACA 63012A AIRFOIL.dat` alpha=-4 nCrit=9 — transition_front_oscillation
- `NACA 63012A AIRFOIL.dat` alpha=+4 nCrit=9 — transition_front_oscillation
- `NACA 8-H-12 AIRFOIL.dat` alpha=-2 nCrit=9 — transition_front_oscillation (rescued)
- `RAE 103 AIRFOIL.dat` alpha=-4 nCrit=9 — transition_front_oscillation (rescued)
- `RAE 103 AIRFOIL.dat` alpha=+6 nCrit=5 — transition_front_oscillation
- `RAE 103 AIRFOIL.dat` alpha=+6 nCrit=9 — transition_front_oscillation (rescued)
- `S2027.dat` alpha=+8 nCrit=9 — transition_front_oscillation (rescued)
- `S4083 (8%).dat` alpha=+6 nCrit=5 — transition_front_oscillation (rescued)
- `SD6080 (9.2%).dat` alpha=+6 nCrit=5 — transition_front_oscillation
- `SD6080 (9.2%).dat` alpha=+8 nCrit=9 — transition_front_oscillation (rescued)
- `SPICA  11.73% smoothed.dat` alpha=+8 nCrit=5 — transition_front_oscillation (rescued)
- `USA 5 AIRFOIL.dat` alpha=-2 nCrit=5 — transition_front_oscillation (rescued)
- `WORTMANN FX 63-137 AIRFOIL.dat` alpha=+6 nCrit=5 — transition_front_oscillation

---

# Analysis (manual sections — `--report` regenerates only the tables above)

## Comparison to the prior 25-foil audit (commit 29e568d, bench/results/REPORT.md)

Same failure taxonomy, no new modes: `transition_front_oscillation` (6.4%) >
`nan_lock` (3.8%) > `no_convergence` (3.2%) > `diverged` (2.4%). Headline
cold-convergence 84.2% vs the prior 85.5–86.1% — consistent given this grid is
all Re=2e6 (the prior grid's Re=2e5/1e6 cells converged at ~90%) and extends to
alpha=+8 (69.2% here; prior alpha=+10 was 70%). The alpha trend reproduces:
best near 0–2 deg (92%), degrading toward high incidence. nCrit is flat
(84.0% vs 84.3%), as before (88/88/86 at 5/9/12).

Differences worth noting:
- **Zero `acoustic_nan` cases** (prior: 1.4%). That mode was a low-Re WPS/Amiet
  degeneracy; this all-Re=2e6 grid never triggers it. The mode is still wired
  and reported separately — it simply doesn't occur here.
- **`newton_iterations`=100 on failures is a sentinel** (unset
  `glob.convergenceIteration`); the solve_coupled cap is 60 iterations and the
  slowest converged case took 59.
- **No crashes**: all 840 cases ran; the per-foil subprocess isolation was never
  needed.

## Light diagnosis of 5 failing cases

`GFOIL_DEBUG=1` only instruments the acoustic chain (`sound.hpp`,
`WPSmodels.hpp`) — there is no residual/omega/transition trace in the aero
solver — so diagnosis used the C++ classifier plus non-invasive probes:
neighbouring alphas (±0.5/±1 deg), rtol=1e-5 retry, warm-start from the nearest
converging neighbour, and geometry statistics.

| case | mode | probe evidence | classification |
|---|---|---|---|
| OAF128, alpha=0, nCrit=5 | nan_lock | **all 14 grid cases fail nan_lock** (both nCrit, every alpha); ±0.5/±1 also fail; rtol-insensitive; geometry smooth (12.6% thick, TE gap 1e-5) | nan_lock, **new presentation: foil-wide/alpha-independent** (prior audit saw isolated points). Geometry-systematic Jacobian breakdown, unrescuable by backstepping (no converging neighbour exists) |
| EPPLER 864 STRUT, alpha=0, nCrit=9 | no_convergence | 38.4% thick strut; nc=9 fails 5/7 alphas with mixed modes (no_conv/tfo/nan_lock), nc=5 mostly converges; rtol-insensitive | thick-section separation regime, out of the integral-BL envelope; matches the prior "genuine BL attractor" class |
| E186 (10.27%), alpha=-2, nCrit=5 | transition_front_oscillation | neighbours converge cold (C18–C41); **warm-start from -1.5 also fails**; rtol-insensitive | genuine multi-node attractor at an isolated point — same class as the documented NACA 0008-34 alpha=-2.6 permanent case |
| EPPLER 399, alpha=+6, nCrit=5 | transition_front_oscillation | neighbours converge (5.0:C27 … 6.5:C43); warm-start "succeeds" only via the it=0 stale accept (see below), i.e. not really | genuine limit cycle at an isolated point; `rescue_conv=0` is correct |
| GOE 117, alpha=+6, nCrit=9 | diverged | failing band 5.0–6.5 (diverged/tfo/nan_lock), 7.0 converges C9; warm-start from 7.0 is an it=0 stale accept | transition-sensitive failure band, attractor class |

## Probable defect found: warm-restart stale accept (driver-level, NOT solver)

While probing warm starts, `run_forward(inp, restart)` was caught returning
`conv=1, newton_iterations=0` with the **donor's solution unchanged** after a
0.5 deg alpha change. Ground-truth case (EPPLER 399, nCrit=5):

```
cold 5.0:            it=27  CL=1.27228   <- truth
cold 5.5:            it=30  CL=1.32060   <- donor
warm 5.0 <- 5.5:     it=0   CL=1.32056   <- donor state accepted as "alpha=5.0"
warm 6.0 <- 6.5:     it=0   CL=1.40844   <- same defect (cold 6.0 doesn't even converge)
warm 5.5 <- 6.5:     it=32  CL=1.32060   <- 1.0 deg jump re-solves correctly
warm 5.0 <- 6.5:     it=33  CL=1.27228   <- 1.5 deg jump re-solves correctly
warm 6.5 <- 6.5:     it=8   CL=1.40848   <- same-alpha restart does NOT instant-accept
```

The stale accept was observed exactly at 0.5 deg jumps — the step size
`standard_run`'s forward-stepping uses — while >=1 deg jumps re-solve. CL error
when it triggers: +0.048 (~3.7%). Impact: any result delivered through the
fwd_run backstepping/forward-stepping path may be up to 0.5 deg stale; the 97
"rescued" cases in this sweep and any production run that printed
"Starting backstepping" are suspect. The same-alpha it=8 vs shifted-alpha it=0
asymmetry suggests the entry-residual convergence check in solve_coupled
interacts with state reindexing (stagpoint_move) in a way that depends on
whether the stagnation node moved — root-causing this is solver work, out of
scope for this measurement pass.

**Proposed safe change (described only, per scope):** in `GFoil/gfoil.py`
(`_call_forward`/`standard_run`), treat a warm-started result with
`newton_iterations == 0` and a changed alpha as NOT converged (fall back to a
cold solve at that alpha, or shrink the step). Python-only, no solver/AD-path
changes, directly evidenced by the table above. A C++ root-cause pass on the
entry residual check (and on `resid_rms`'s Rsize=3*(Ncoords+Nwake) row
coverage) should follow separately before trusting any backstepped result.

## Baseline contract

Future solver changes (triplet scatter-add, factorization swaps, etc.) must
reproduce, on `python3 bench/cold_start_sweep.py` (venv
/home/pa20830/NewGradientVal/venv, seed=0 sample):

- cold-converged >= 707/840, with no new failure modes beyond the four above;
- the per-alpha profile within ~1pt (in particular alpha=+8 ~ 69%);
- iteration distribution: median 16, p90 41, max < 60 on converged cases;
- the full failing-case list (above) as the expected-failure set — any newly
  failing (foil, alpha, nCrit) is a regression even if the totals match.

---

# Addendum (Part A): rescue re-measurement with warm-restart guard

`GFoil/gfoil.py` now rejects a warm-restart iteration-0 accept at a changed alpha (`_is_stale_warm_accept`) and cold-retries the same alpha. Re-running the full fwd_run continuation on the 133 cold failures (`bench/rescue_guard.py`, guard active):

- prior rescues (unguarded `rescue_conv`): **97/133**
- corrected rescues (`rescue_conv_guarded`): **87/133**
- prior "rescues" that were stale accepts (lost): **13**
- newly-rescued under the guard: **3**
- guard fired during **32** of the 133 rescue runs

The corrected count replaces the suspect 97. `rescue_conv` is kept as-is in the CSV; `rescue_conv_guarded` is the new column.

## Cases that flipped 1 -> 0 (prior rescue was a stale accept)

- `AH 79-100 C AIRFOIL.dat` alpha=-2 nCrit=9 — transition_front_oscillation
- `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` alpha=+0 nCrit=5 — transition_front_oscillation
- `EPPLER 1098 AIRFOIL.dat` alpha=-4 nCrit=5 — transition_front_oscillation
- `EPPLER 502 AIRFOIL.dat` alpha=-4 nCrit=9 — nan_lock
- `FX 61-168 AIRFOIL.dat` alpha=-4 nCrit=5 — no_convergence
- `FX 61-168 AIRFOIL.dat` alpha=-4 nCrit=9 — nan_lock
- `FX60-100 10.0% smoothed.dat` alpha=-2 nCrit=5 — nan_lock
- `FX60-100 10.0% smoothed.dat` alpha=+4 nCrit=5 — transition_front_oscillation
- `GOE 425 AIRFOIL.dat` alpha=+0 nCrit=5 — transition_front_oscillation
- `GOE 513 AIRFOIL.dat` alpha=-4 nCrit=5 — no_convergence
- `GOE 770 AIRFOIL.dat` alpha=-4 nCrit=5 — transition_front_oscillation
- `SD6080 (9.2%).dat` alpha=+0 nCrit=5 — nan_lock
- `USA 5 AIRFOIL.dat` alpha=-2 nCrit=5 — transition_front_oscillation

## Cases that flipped 0 -> 1 under the guard

- `GOE 117 (MVA MK.4) AIRFOIL.dat` alpha=+6 nCrit=9 — diverged
- `GOE 328 AIRFOIL.dat` alpha=+2 nCrit=5 — transition_front_oscillation
- `RAE 103 AIRFOIL.dat` alpha=+6 nCrit=5 — transition_front_oscillation

## Full per-case rescue table (cold failures)

| foil | alpha | nCrit | mode | rescue_conv | rescue_conv_guarded | guard_fired |
|---|---|---|---|---|---|---|
| `AH 79-100 C AIRFOIL.dat` | +8 | 9 | diverged | 1 | 1 | 0 |
| `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` | +6 | 9 | diverged | 1 | 1 | 0 |
| `CH10 (smoothed).dat` | +8 | 5 | diverged | 1 | 1 | 0 |
| `EPPLER 1098 AIRFOIL.dat` | +0 | 5 | diverged | 1 | 1 | 0 |
| `FX 61-168 AIRFOIL.dat` | +8 | 5 | diverged | 1 | 1 | 0 |
| `FX60-100 10.0% smoothed.dat` | +8 | 9 | diverged | 1 | 1 | 1 |
| `GOE 117 (MVA MK.4) AIRFOIL.dat` | +6 | 9 | diverged | 0 | 1 | 1 |
| `GOE 184 (MVA H.29) AIRFOIL.dat` | -4 | 9 | diverged | 1 | 1 | 0 |
| `GOE 346 (FRIEDRICHSHAFEN-STAAKEN) AIRFOIL.dat` | +8 | 5 | diverged | 1 | 1 | 0 |
| `GOE 393 AIRFOIL.dat` | -4 | 9 | diverged | 1 | 1 | 0 |
| `GOE 458 AIRFOIL.dat` | +4 | 9 | diverged | 1 | 1 | 0 |
| `GOE 513 AIRFOIL.dat` | +8 | 9 | diverged | 1 | 1 | 0 |
| `GOE 596 AIRFOIL.dat` | +8 | 9 | diverged | 1 | 1 | 0 |
| `GOE 777 AIRFOIL.dat` | +8 | 9 | diverged | 1 | 1 | 0 |
| `LWK 80-150-K25.dat` | +4 | 5 | diverged | 0 | 0 | 1 |
| `MH 121  8.76%.dat` | -2 | 5 | diverged | 1 | 1 | 0 |
| `MH 121  8.76%.dat` | +2 | 9 | diverged | 1 | 1 | 0 |
| `NASA SC(2)-0012 AIRFOIL.dat` | -4 | 5 | diverged | 1 | 1 | 0 |
| `NASA SC(2)-0012 AIRFOIL.dat` | +4 | 5 | diverged | 1 | 1 | 0 |
| `S9033 (7.5%).dat` | +8 | 5 | diverged | 1 | 1 | 0 |
| `AH 79-100 C AIRFOIL.dat` | +6 | 5 | nan_lock | 1 | 1 | 0 |
| `CH10 (smoothed).dat` | +8 | 9 | nan_lock | 1 | 1 | 0 |
| `E374.dat` | +4 | 9 | nan_lock | 1 | 1 | 0 |
| `E374.dat` | +8 | 5 | nan_lock | 1 | 1 | 0 |
| `EPPLER 502 AIRFOIL.dat` | -4 | 9 | nan_lock | 1 | 0 | 1 |
| `EPPLER 502 AIRFOIL.dat` | +8 | 5 | nan_lock | 1 | 1 | 0 |
| `EPPLER 864 STRUT AIRFOIL.dat` | -2 | 5 | nan_lock | 1 | 1 | 0 |
| `EPPLER 864 STRUT AIRFOIL.dat` | +6 | 9 | nan_lock | 0 | 0 | 0 |
| `FX 61-168 AIRFOIL.dat` | -4 | 9 | nan_lock | 1 | 0 | 1 |
| `FX60-100 10.0% smoothed.dat` | -2 | 5 | nan_lock | 1 | 0 | 1 |
| `GOE 443 AIRFOIL.dat` | -4 | 5 | nan_lock | 1 | 1 | 0 |
| `GOE 443 AIRFOIL.dat` | +4 | 5 | nan_lock | 1 | 1 | 0 |
| `GOE 513 AIRFOIL.dat` | +6 | 9 | nan_lock | 1 | 1 | 0 |
| `GOE 777 AIRFOIL.dat` | +8 | 5 | nan_lock | 1 | 1 | 0 |
| `MH 121  8.76%.dat` | +2 | 5 | nan_lock | 0 | 0 | 0 |
| `NACA 2415.dat` | +8 | 5 | nan_lock | 1 | 1 | 0 |
| `OAF128 AIRFOIL.dat` | -4 | 5 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | -4 | 9 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | -2 | 5 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | -2 | 9 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +0 | 5 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +0 | 9 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +2 | 5 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +2 | 9 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +4 | 5 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +4 | 9 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +6 | 5 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +6 | 9 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +8 | 5 | nan_lock | 0 | 0 | 0 |
| `OAF128 AIRFOIL.dat` | +8 | 9 | nan_lock | 0 | 0 | 0 |
| `SD6080 (9.2%).dat` | +0 | 5 | nan_lock | 1 | 0 | 1 |
| `USA 5 AIRFOIL.dat` | +6 | 9 | nan_lock | 1 | 1 | 0 |
| `AH 79-100 C AIRFOIL.dat` | +0 | 9 | no_convergence | 1 | 1 | 0 |
| `CH10 (smoothed).dat` | -2 | 9 | no_convergence | 1 | 1 | 1 |
| `EPPLER 1098 AIRFOIL.dat` | +8 | 5 | no_convergence | 1 | 1 | 0 |
| `EPPLER 502 AIRFOIL.dat` | +6 | 9 | no_convergence | 1 | 1 | 0 |
| `EPPLER 582 AIRFOIL.dat` | +8 | 9 | no_convergence | 1 | 1 | 0 |
| `EPPLER 864 STRUT AIRFOIL.dat` | -2 | 9 | no_convergence | 0 | 0 | 0 |
| `EPPLER 864 STRUT AIRFOIL.dat` | +0 | 9 | no_convergence | 0 | 0 | 0 |
| `EPPLER 864 STRUT AIRFOIL.dat` | +4 | 5 | no_convergence | 1 | 1 | 0 |
| `FX 61-168 AIRFOIL.dat` | -4 | 5 | no_convergence | 1 | 0 | 1 |
| `FX 67-K-150-17 AIRFOIL.dat` | +6 | 5 | no_convergence | 1 | 1 | 0 |
| `FX 67-K-150-17 AIRFOIL.dat` | +8 | 9 | no_convergence | 0 | 0 | 1 |
| `FX60-100 10.0% smoothed.dat` | +6 | 9 | no_convergence | 1 | 1 | 0 |
| `GOE 226 (MVA H.36) AIRFOIL.dat` | -2 | 9 | no_convergence | 1 | 1 | 1 |
| `GOE 513 AIRFOIL.dat` | -4 | 5 | no_convergence | 1 | 0 | 1 |
| `GOE 513 AIRFOIL.dat` | -4 | 9 | no_convergence | 1 | 1 | 0 |
| `GOE 602 MOD. AIRFOIL.dat` | +6 | 9 | no_convergence | 1 | 1 | 0 |
| `GOE 629 AIRFOIL.dat` | -2 | 5 | no_convergence | 1 | 1 | 0 |
| `GOE 629 AIRFOIL.dat` | +8 | 9 | no_convergence | 1 | 1 | 1 |
| `GOE 770 AIRFOIL.dat` | +6 | 9 | no_convergence | 1 | 1 | 0 |
| `GOE 777 AIRFOIL.dat` | +0 | 9 | no_convergence | 1 | 1 | 0 |
| `HQ 2.5-8 AIRFOIL.dat` | +6 | 5 | no_convergence | 1 | 1 | 0 |
| `LWK 80-150-K25.dat` | -4 | 5 | no_convergence | 0 | 0 | 1 |
| `MH 121  8.76%.dat` | -2 | 9 | no_convergence | 1 | 1 | 0 |
| `NACA 2415.dat` | -4 | 9 | no_convergence | 1 | 1 | 0 |
| `NACA 65-210.dat` | +8 | 9 | no_convergence | 1 | 1 | 0 |
| `S4083 (8%).dat` | +8 | 5 | no_convergence | 1 | 1 | 0 |
| `WORTMANN FX 63-137 AIRFOIL.dat` | +8 | 9 | no_convergence | 1 | 1 | 0 |
| `AG17.dat` | +6 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `AG17.dat` | +6 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `AH 79-100 C AIRFOIL.dat` | -2 | 9 | transition_front_oscillation | 1 | 0 | 1 |
| `AH21 7% version (Andrew Hollom).dat` | -4 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` | +0 | 5 | transition_front_oscillation | 1 | 0 | 1 |
| `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` | +2 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `BELL-WORTMANN FX 69-H-083 AIRFOIL.dat` | +6 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `E186 (10.27%).dat` | -2 | 5 | transition_front_oscillation | 0 | 0 | 0 |
| `EPPLER 1098 AIRFOIL.dat` | -4 | 5 | transition_front_oscillation | 1 | 0 | 1 |
| `EPPLER 399 AIRFOIL.dat` | +6 | 5 | transition_front_oscillation | 0 | 0 | 0 |
| `EPPLER 399 AIRFOIL.dat` | +8 | 5 | transition_front_oscillation | 0 | 0 | 0 |
| `EPPLER 399 AIRFOIL.dat` | +8 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `EPPLER 502 AIRFOIL.dat` | +4 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `EPPLER 502 AIRFOIL.dat` | +8 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `EPPLER 637 AIRFOIL.dat` | +0 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `EPPLER 864 STRUT AIRFOIL.dat` | +2 | 9 | transition_front_oscillation | 0 | 0 | 0 |
| `EPPLER 864 STRUT AIRFOIL.dat` | +4 | 9 | transition_front_oscillation | 0 | 0 | 0 |
| `FX 61-168 AIRFOIL.dat` | +8 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `FX 67-K-150-17 AIRFOIL.dat` | +8 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `FX60-100 10.0% smoothed.dat` | +4 | 5 | transition_front_oscillation | 1 | 0 | 1 |
| `FX60-100 10.0% smoothed.dat` | +4 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `FX60-100 10.0% smoothed.dat` | +8 | 5 | transition_front_oscillation | 1 | 1 | 1 |
| `GOE 117 (MVA MK.4) AIRFOIL.dat` | +6 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `GOE 117 (MVA MK.4) AIRFOIL.dat` | +8 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `GOE 328 AIRFOIL.dat` | +2 | 5 | transition_front_oscillation | 0 | 1 | 1 |
| `GOE 346 (FRIEDRICHSHAFEN-STAAKEN) AIRFOIL.dat` | +8 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `GOE 393 AIRFOIL.dat` | -4 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `GOE 425 AIRFOIL.dat` | +0 | 5 | transition_front_oscillation | 1 | 0 | 1 |
| `GOE 443 AIRFOIL.dat` | +6 | 5 | transition_front_oscillation | 0 | 0 | 1 |
| `GOE 443 AIRFOIL.dat` | +8 | 5 | transition_front_oscillation | 0 | 0 | 1 |
| `GOE 770 AIRFOIL.dat` | -4 | 5 | transition_front_oscillation | 1 | 0 | 1 |
| `GOE 770 AIRFOIL.dat` | -2 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `GOE 777 AIRFOIL.dat` | -4 | 9 | transition_front_oscillation | 1 | 1 | 1 |
| `GOE 777 AIRFOIL.dat` | -2 | 5 | transition_front_oscillation | 1 | 1 | 1 |
| `HQ 1.5-8 AIRFOIL.dat` | +2 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `HQ 1.5-8 AIRFOIL.dat` | +2 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `HQ 1.5-8 AIRFOIL.dat` | +4 | 5 | transition_front_oscillation | 0 | 0 | 0 |
| `HQ 1.5-8 AIRFOIL.dat` | +6 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `HQ 2.5-8 AIRFOIL.dat` | -2 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `LWK 80-150-K25.dat` | +8 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `MH 49.dat` | +8 | 9 | transition_front_oscillation | 1 | 1 | 1 |
| `NACA 63012A AIRFOIL.dat` | -4 | 9 | transition_front_oscillation | 0 | 0 | 1 |
| `NACA 63012A AIRFOIL.dat` | +4 | 9 | transition_front_oscillation | 0 | 0 | 1 |
| `NACA 8-H-12 AIRFOIL.dat` | -2 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `RAE 103 AIRFOIL.dat` | -4 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `RAE 103 AIRFOIL.dat` | +6 | 5 | transition_front_oscillation | 0 | 1 | 1 |
| `RAE 103 AIRFOIL.dat` | +6 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `S2027.dat` | +8 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `S4083 (8%).dat` | +6 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `SD6080 (9.2%).dat` | +6 | 5 | transition_front_oscillation | 0 | 0 | 1 |
| `SD6080 (9.2%).dat` | +8 | 9 | transition_front_oscillation | 1 | 1 | 0 |
| `SPICA  11.73% smoothed.dat` | +8 | 5 | transition_front_oscillation | 1 | 1 | 0 |
| `USA 5 AIRFOIL.dat` | -2 | 5 | transition_front_oscillation | 1 | 0 | 1 |
| `WORTMANN FX 63-137 AIRFOIL.dat` | +6 | 5 | transition_front_oscillation | 0 | 0 | 0 |


---

# Addendum (Part B): solve_coupled entry-accept root cause + fix

## Phase 1 — instrumented confirmation (GFOIL_DEBUG, `bench/b1_instrument.py`)

EPPLER 399, nCrit=5, Re=2e6, rtol=1e-6. `solve_coupled` prints BL-row RMS
(existing `residualNorm`), ue-row RMS (filled via `ue_residual_kernel` into a
scratch — write-not-accumulate, refilled identically by the later `solve_glob`),
`stagIndex`, and `Is` sizes at iteration 0 (ENTRY) and at convergence.

| case | entry BL_rms | entry ue_rms | rtol | it | outcome |
|---|---|---|---|---|---|
| (a) warm 5.0 <- 5.5 | **1.80e-07** | **0.0422** | 1e-6 | 0 | instant accept, donor CL echoed |
| (b) warm 6.5 <- 6.5 | 0.621 | 0.0098 | 1e-6 | 8 | re-settles (see anomaly) |
| (c) cold 5.0 (control) | 0.0374 | 0.0962 | 1e-6 | 27 | both large, normal solve |

**Smoking gun (case a) — predicted signature confirmed.** At a warm entry the BL
rows are the donor's converged residuals (1.8e-7 << rtol) while the ue-coupling
rows — where the alpha change enters (via `uewi`/`gammas` in
`ue_residual_kernel`, evaluated only inside `solve_glob` AFTER the test) — carry
a residual of 0.042, four orders of magnitude above rtol. `resid_rms(glob.R,
Rsize)` with `Rsize = 3*(Ncoords+Nwake)` covers only the BL rows, so the alpha
mismatch is invisible and the donor state is accepted at iteration 0.

## Same-alpha it=8 anomaly (case b) — a SECOND, distinct defect (scoped, not fixed)

The donor (cold 6.5) **converged** at `stagIndex=[82,83]` (Is sizes 83/117). The
warm 6.5<-6.5 re-entry **enters** at `stagIndex=[81,82]` (Is sizes 82/118), with
`xi` shifted (0.9881 vs the donor's converged 0.9922) and entry `BL_rms=0.62`.
Root cause: on a warm restart `runCode` runs the INVISCID `stagpoint_find`
(+`identify_surfaces`/`set_wake_gap`), loads the donor states, then a SINGLE
viscous `stagpoint_move` — which lands one node off the donor's converged
stagnation index. That reindexes the BL stations and shifts `distFromStag`, so
the BL rows are no longer the donor's converged residuals (hence no instant
accept; 8 iterations to walk back to [82,83]). `RestartState.stag` is stored but
never re-imposed on warm entry — `stagpoint_move` always recomputes it. The
0.5-deg vs same-alpha asymmetry in the original table is exactly this: when the
inviscid-seeded `stagpoint_move` happens to land on the donor's converged node
(case a) the BL rows read as converged and it=0; when it lands off-by-one
(case b) they don't. wgap was identical across donor/warm (same geometry+alpha →
same inviscid stag), so it is NOT the differing quantity. **Scope:** documented,
not fixed in this pass; a fix would re-impose `RestartState.stag` (and refresh
`set_wake_gap`) on warm entry.

## Phase 2 — two fixes behind `GFOIL_BFIX` (default Fix B)

- **Fix A** (`GFOIL_BFIX=A`): honest criterion `max(BL_rms, ue_rms) < rtol`
  (the ue rows filled before the test; a max, not a pooled 4*Nsys RMS, to keep
  the per-equation rtol semantics of the BL rows).
- **Fix B** (default / `GFOIL_BFIX=B`): skip the iteration-0 accept on a warm
  (restart) entry, forcing >=1 Newton iteration (`warmEntry` threaded from
  `runCode`). Cold paths bit-identical by construction.
- `GFOIL_BFIX=off`: original BL-only criterion (defect reproduction).

### Measurements

**EPPLER 399 warm table** (`bench/b3_eppler.py`, direct `_call_forward`):

| mode | warm 5.0 <- 5.5 | CL | verdict |
|---|---|---|---|
| off | it=0 | 1.32056 | STALE (defect) |
| A | it=26 | 1.27228 | cold-truth |
| B | it=26 | 1.27228 | cold-truth |

**Golden regression** (`tests/regression_test.py --test`): Fix B **22/22
bit-identical** (max rel_err 0). Fix A **also 22/22 bit-identical** (rel_err 0) —
on the golden case the ue rows are already converged when the BL rows converge,
so the max-criterion crosses rtol at the same iteration. Goldens NOT regenerated.

**840-case cold sweep** (`bench/cold_only_sweep.py`, seed=0 sample) vs the
original baseline (`cold_baseline_orig.csv`):

| metric | baseline | Fix B | Fix A |
|---|---|---|---|
| cold-converged | 707/840 | 707/840 | 707/840 |
| conv-flag diffs | — | 0 | 0 |
| failure_mode diffs | — | 0 | 0 |
| newton_iter diffs | — | 0 | 0 |
| over 60-iter cap | 0 | 0 | 0 |
| iters median/p90/max | 16/41/59 | 16/41/59 | 16/41/59 |
| per-alpha (+8) | 83/120 | 83/120 | 83/120 |

Both fixes leave **every** cold trajectory bit-identical to the baseline. No
newly-failing case under Fix A ⇒ the B.3(iii) spot-check (re-run newly-failing
cases at rtol=1e-8 under the old criterion) has no candidates: on this cold grid
the old criterion never declared convergence on an unconverged ue system.

**Wall-time** (30-rep forward median, golden case): off ≈90 ms, Fix B ≈90 ms
(no-op on cold), Fix A ≈92 ms (~2 %, the extra per-iteration ue-kernel).

### Decision

**Fix B active by default; Fix A retained behind `GFOIL_BFIX=A`.** Fix B is
sufficient for the observed defect, provably bit-identical on cold paths, and
zero-cost. Fix A is the more principled criterion and is empirically
bit-identical on the golden + 840 grid, but universal bit-identity is not
guaranteed, so it stays behind the switch (promotable as a one-liner). **Baseline
contract unchanged** — Part C gates against `cold_sweep_B.csv` (== original
baseline).

---

# Addendum (Part C): scatter-add Newton linear solve (`solve_sys_sparse`)

Replaces the per-iteration `setFromTriplets` rebuild with a pattern-keyed
scatter-add into a cached, compressed matrix. Keyed to the SAME FNV pattern hash
already used for the analyzePattern-once cache: on a pattern change, build A with
`setFromTriplets`, `analyzePattern`, and record `slot[k]` = the `valuePtr()` index
for triplet k (binary search of the row within its column's compressed inner
range) plus `is_first[k]` (first triplet to touch each slot). On an unchanged
pattern, skip `setFromTriplets`: the first triplet per slot ASSIGNS (seeds
verbatim), duplicates ADD in ascending-k order.

## C.1 pattern churn (`bench/c1_pattern.py`, fresh subprocess per case)

| case | solves | pattern changes | scatter (fast path) |
|---|---|---|---|
| golden NACA0012 a2 nc5 | 9 | 1 (11%) | 8 (89%) |
| S4083 (8%) +6 nc9 | 59 | 8 (14%) | 51 (86%) |
| GOE 328 +8 nc9 | 58 | 5 (9%) | 53 (91%) |
| GOE 303 +0 nc9 | 58 | 12 (21%) | 46 (79%) |
| WORTMANN FX 63-137 +6 nc9 | 57 | 12 (21%) | 45 (79%) |
| AH 79-100 C +6 nc9 | 56 | 18 (32%) | 38 (68%) |

Even the worst slow case keeps 68% of iterations on the scatter path; none exceed
~half churn. The win is preserved once stag/transition settle.

## C.2 bit-identity (mandatory, GFOIL_CVERIFY — every iteration)

Build A both ways and memcmp `outerIndexPtr`/`innerIndexPtr`/`valuePtr` each call.
Initial mismatch was **signed zero**: the triplet list is duplicate-free here
(`nonZeros() == nnz` in every case — `addColumnValues` pre-sums in place), so
`setFromTriplets` stores each value verbatim including `-0.0`, whereas a
`fill(0.0)+=` flips `-0.0` to `+0.0` (`0.0 + -0.0 == +0.0`). Fixed by ASSIGN-first
(seed verbatim) / ADD-duplicates. After the fix: **0 memcmp failures across all
297 solves** of the golden + 5 slow cases, including the golden AD recording pass.
`dup_iters = 0` everywhere (no case has duplicate (row,col) triplets), and the
`setFromTriplets` determinism self-check (`ref_det=1`) held throughout. The EF's
`data->A` is the same cached matrix (C.3), so the adjoint sees identical values.

## C.4 gates

- (i)+(ii) Full regression (`tests/regression_test.py --test`): **22/22
  bit-identical** with scatter (default) and with `GFOIL_NOSCATTER=1` (control).
- (iii) 840-case cold sweep vs the post-B baseline (`cold_sweep_B.csv`):
  **0 conv-flag, 0 failure_mode, 0 newton_iteration diffs**; 707/840, median 16,
  p90 41, max 59 — identical Newton trajectories.
- (iv) Wall-time (`bench/c4_timing.py`, scatter vs `GFOIL_NOSCATTER`, same binary,
  full `run_forward`):

| case | noscatter | scatter | delta |
|---|---|---|---|
| golden NACA0012 a2 nc5 | 89.7 ms | 86.2 ms | -3.9% |
| S4083 (8%) +6 nc9 | 448.9 | 410.6 | -8.5% |
| BELL-WORTMANN FX69 +4 nc5 | 430.5 | 388.7 | -9.7% |
| GOE 303 +0 nc9 | 472.0 | 461.6 | -2.2% |
| GOE 328 +8 nc9 | 468.6 | 421.6 | -10.0% |
| GOE 458 -2 nc9 | 461.6 | 402.4 | -12.8% |
| NASA MS(1)-0313 +8 nc9 | 443.5 | 415.0 | -6.4% |
| WORTMANN FX 63-137 +6 nc9 | 447.1 | 418.4 | -6.4% |
| AH 79-100 C +6 nc9 | 443.8 | 429.5 | -3.2% |
| MH 93 15.98% -4 nc5 | 435.0 | 405.9 | -6.7% |
| E374 -2 nc9 | 443.9 | 414.4 | -6.6% |

High-iteration stable-pattern cases gain most (GOE 458 -12.8%, GOE 328 -10.0%);
the highest-churn case (AH 79-100 C, 32% rebuilds) gains least (-3.2%), and the
golden (few iterations, acoustic-dominated total) -3.9%. Consistent with C.1.

---

# Addendum (Part B second defect, fixed): warm-restart stagnation reindex

The same-alpha it=8 anomaly scoped in the Part B addendum is now fixed. On a warm
entry `runCode` ran `identify_surfaces`/`set_wake_gap`/`calc_ue_m` from the
INVISCID `stagpoint_find`, then a single viscous `stagpoint_move` whose sign-scan
is seeded from `isol.stagIndex`. Left at the inviscid value it landed one node off
the donor's converged stag, reindexing the BL stations. Fix: seed the bracket from
the donor's converged stag (`RestartState.stag`, now plumbed through the pybind
warm path and read at load); `stagpoint_move`'s `identify_surfaces` then rebuilds
`Is`/`distFromStag` to the donor configuration. `wgap` is already donor-identical
(keyed to the inviscid stag, same for the same geometry+alpha). `GFOIL_NOSTAGSEED`
keeps the pre-fix behaviour for A/B.

## Phase 1 — mechanism (GFOIL_DEBUG, `bench/stagseed_probe.py`)

EPPLER 399 nCrit=5, warm 6.5<-6.5 (donor cold-6.5 converges at stag [82,83]):

| mode | inviscid stag | post-move stag | entry BL_rms | Is sizes | it |
|---|---|---|---|---|---|
| NOSTAGSEED (pre-fix) | [81,82] | [81,82] | 0.621 | 82/118 | 8 |
| seed (fix) | [81,82] | **[82,83]** | **7.7e-09** | **83/117** | **1** |

The seeded post-move stag, Is sizes (83/117) and `distFromStag` (xi=0.992173) match
the donor's converged values exactly — it reproduces the donor configuration, not
a new one. (a) The discrepancy is purely the sign-scan seed; the loaded states are
consistent with the donor stag, so seeding there yields BL_rms~1e-8. (b)
`stagpoint_move` rebuilds `Is` every iteration via `identify_surfaces`, so the
donor's converged Is(83/117) tracks its converged stag [82,83]; the seed
reproduces it. For warm 5.0<-5.5 the inviscid stag already equals the donor
[83,84], so the seed is a no-op there.

## Gates

1. Regression: **22/22 bit-identical** (golden is cold; the warm seed path is not
   taken — cold runs never have a warmStart).
2. EPPLER same-alpha probe: warm 6.5<-6.5 **8 -> 1 iteration**; warm 5.0<-5.5 stays
   a real re-solve under Fix B (**it=26, nonzero**) with cold-truth **CL 1.27228**
   (|Δ|=2.8e-6) — the stale-accept fix remains effective.
3. 133-case rescue (`bench/rescue_cost.py`, seed vs NOSTAGSEED on the same binary,
   both with Fix B + scatter): rescue rate **94/133 -> 94/133** (unchanged, no
   regressions); converged-iteration cost over the 94 commonly-rescued cases
   **5157 -> 4408 (-14.5%)**; 70 cases cheaper, 8 marginally costlier. Largest
   wins: AH 79-100 C +8/nc9 142->87, GOE 346 +8/nc5 59->19, GOE 117 +6/nc5 89->54.
4. 840 cold sweep vs the post-C baseline: **0 conv / 0 failure_mode / 0
   newton_iteration diffs** — warm path untouched on cold runs.
