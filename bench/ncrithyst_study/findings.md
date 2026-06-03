# ncrithyst convergence study — findings

*Seed=0, 80 aerofoils from `Smoothed_TEfixed_linear/`. Grid: ncrithyst[0.0, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0] x Re{2e5,5e5,1e6,3e6,6e6,9e6} x alpha{0,2,4,6,8,10,12}° x nCrit{6,9,12}, FREE transition only. Bare single Newton solve (`_build_input_dict`+`cpp.run_forward`, no backstepping/repanel) so the continuation fallback cannot mask the effect. Total runs: 90,720. Built from a HEAD worktree that still exposes ncrithyst; no library source modified.*

## Headline

**Qualified yes — but it is not a usable convergence knob.** Raising `ncrithyst` raises the bare-Newton convergence rate monotonically (78.0% at 0 → 90.5% at 3.0), but the gain is bought by **over-damping the transition front**: beyond ~0.8 it changes the converged solution (systematic drag rise) rather than stabilising it. At the physically safe small values (≤0.4, incl. the historical default **0.2**), the effect on convergence is **within noise** (78.05% → 78.38%). So there is no value that materially improves convergence *without* perturbing the answer. Recommendation below: do **not** inflate it.

## Convergence rate vs ncrithyst (full curve)

| ncrithyst | 0 | 0.1 | 0.2 | 0.4 | 0.8 | 1.2 | 1.6 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|
| aero conv % | 78.05 | 78.45 | 78.38 | 78.02 | 84.28 | 87.09 | 88.41 | 89.14 | 90.47 |

Flat-then-rising: 0.1–0.4 are indistinguishable from OFF; the gain only appears once the margin is wide enough to catch the node (≥0.8) and then keeps climbing with no plateau inside [0, 3.0].

## The catch — solution drift on matched-converged cases

On the 6,030 cases that converge at EVERY ncrithyst, the converged solution drifts systematically (ΔCD is consistently **positive** — front pinned forward → more turbulent → higher drag):

| ncrithyst | aero conv % | mean \|ΔCL\| | mean ΔCD (signed) | % cases CL moved >1e-4 |
|---|---|---|---|---|
| 0 | 78.05 | 0.00e+00 | +0.00e+00 | 0.0 |
| 0.1 | 78.45 | 1.86e-04 | +5.18e-06 | 7.0 |
| 0.2 | 78.38 | 4.00e-04 | +2.77e-06 | 16.7 |
| 0.4 | 78.02 | 1.30e-03 | +2.65e-05 | 37.1 |
| 0.8 | 84.28 | 4.52e-03 | +1.82e-04 | 59.1 |
| 1.2 | 87.09 | 6.13e-03 | +2.85e-04 | 63.8 |
| 1.6 | 88.41 | 6.61e-03 | +3.15e-04 | 64.7 |
| 2 | 89.14 | 6.79e-03 | +3.24e-04 | 65.0 |
| 3 | 90.47 | 6.99e-03 | +3.29e-04 | 65.1 |

A genuine de-oscillation fix would converge to a margin-independent answer (CL/CD plateau once the margin is 'big enough'). Instead CL/CD keep drifting and drag rises monotonically — the hallmark of suppressing *legitimate* front motion, not just spurious 1-node chatter. See `plot_tradeoff.png`.

## Failure modes per ncrithyst

| ncrithyst | diverged | nan_lock | no_convergence | transition_front_oscillation |
|---|---|---|---|---|
| 0 | 428 | 381 | 648 | 756 |
| 0.1 | 433 | 365 | 635 | 739 |
| 0.2 | 435 | 359 | 638 | 747 |
| 0.4 | 456 | 374 | 651 | 735 |
| 0.8 | 228 | 272 | 494 | 591 |
| 1.2 | 174 | 255 | 351 | 521 |
| 1.6 | 132 | 242 | 313 | 481 |
| 2 | 133 | 253 | 269 | 440 |
| 3 | 128 | 249 | 216 | 368 |

`transition_front_oscillation` (the exact mode the hysteresis targets): 756 at h=0 → 747 at h=0.2 (barely moved) → 591 at 0.8 → 368 at 3.0. ALL modes (diverged, nan_lock, no_convergence) fall together at large h, consistent with the front being pinned so the solve simply stops exploring.

## Rescues vs regressions (aero) relative to ncrithyst=0

| ncrithyst | rescued (conv only at h) | regressed (conv only at 0) |
|---|---|---|
| 0.1 | 532 | 491 |
| 0.2 | 726 | 692 |
| 0.4 | 916 | 919 |
| 0.8 | 1299 | 671 |
| 1.2 | 1517 | 605 |
| 1.6 | 1610 | 565 |
| 2 | 1652 | 534 |
| 3 | 1724 | 472 |

Even at large h there is real two-way churn (hundreds of regressions), i.e. for some foils the extra damping *prevents* a previously-good solve.

## Speed (matched-converged Newton iterations)

| ncrithyst | 0 | 0.1 | 0.2 | 0.4 | 0.8 | 1.2 | 1.6 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|
| mean iters | 17.3 | 17.0 | 17.0 | 16.6 | 14.5 | 13.2 | 12.6 | 12.3 | 12.1 |
| median | 12 | 12 | 12 | 13 | 11 | 10 | 9 | 9 | 9 |

## Regime concentration

### aero conv % by Re

| Re | h=0 | h=0.1 | h=0.2 | h=0.4 | h=0.8 | h=1.2 | h=1.6 | h=2 | h=3 |
|---|---|---|---|---|---|---|---|---|---|
| 200000 | 80 | 82 | 81 | 81 | 74 | 71 | 72 | 72 | 75 |
| 500000 | 89 | 89 | 88 | 88 | 88 | 88 | 90 | 92 | 94 |
| 1e+06 | 84 | 84 | 84 | 83 | 86 | 93 | 94 | 95 | 96 |
| 3e+06 | 74 | 75 | 74 | 75 | 87 | 90 | 92 | 92 | 93 |
| 6e+06 | 72 | 71 | 73 | 70 | 86 | 90 | 91 | 92 | 93 |
| 9e+06 | 68 | 69 | 71 | 72 | 84 | 89 | 91 | 91 | 93 |

### aero conv % by alpha

| alpha | h=0 | h=0.1 | h=0.2 | h=0.4 | h=0.8 | h=1.2 | h=1.6 | h=2 | h=3 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 90 | 89 | 89 | 89 | 90 | 89 | 89 | 89 | 89 |
| 2 | 87 | 87 | 89 | 87 | 88 | 88 | 88 | 88 | 89 |
| 4 | 83 | 85 | 84 | 83 | 84 | 86 | 87 | 88 | 89 |
| 6 | 81 | 80 | 81 | 80 | 86 | 88 | 89 | 90 | 91 |
| 8 | 73 | 74 | 74 | 71 | 83 | 87 | 90 | 90 | 92 |
| 10 | 67 | 67 | 67 | 68 | 81 | 86 | 89 | 90 | 92 |
| 12 | 65 | 67 | 66 | 68 | 78 | 84 | 88 | 89 | 91 |

The benefit concentrates at **high Re and high alpha (near stall)** — exactly where transition-front oscillation is expected — but so does the pinning side-effect.

## Recommendation

- **Do not use `ncrithyst` as a convergence fix / do not inflate it.** The large-value convergence gains are over-damping artefacts that bias the converged solution (drag rises monotonically with no plateau).

- **At physically-safe small values (≤0.4) it is essentially inert** for both convergence and solution. The historical default **0.2** is a harmless, near-no-op anti-chatter guard (16.7% of matched cases move CL by >1e-4, mean |ΔCL|≈4e-4) and is a defensible default — but it does **not** measurably help convergence (+0.3 pts, within noise).

- **Real cold-start robustness should come from the continuation / backstepping machinery** (`standard_run`, panel re-distribution), not from widening this hysteresis margin.

- If forced to pick a single hard-coded value on convergence grounds alone, **0.2** is the right conservative choice; anything ≥0.8 trades solution fidelity for a convergence number and should be rejected.
