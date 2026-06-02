# Phase 2 — residual line search: investigated and rejected

## Summary

A residual-monotonicity (Armijo) backtracking line search around the state
update was implemented and benchmarked in two forms. **Both failed to improve
cold-start convergence over Phase 1**, and the strict form badly regressed it.
The line search was therefore **reverted**. The C++ tree after Phase 2 is
identical to Phase 1; this is a documented negative result.

| variant | cold-converged | median iters | wall-time |
|---|---|---|---|
| Phase 1 (no line search) | 1184/1375 = 86.1 % | 13 | 337 s |
| Phase 2 strict Armijo | 1072/1375 = 78.0 % | 12 | 547 s |
| Phase 2 non-monotone (Grippo) | 1183/1375 = 86.0 % | 13 | 547 s |

## What was built

- `update_state.cpp` refactored into `compute_limiter_omega` (the six physical
  limiters → ω) and `apply_state_update(scale)` (the `U += scale·dU` increment
  plus the Hk / negative-ctau / amp-clamp fixes), with a bit-identical
  `update_state` wrapper.
- A line search in `coupled.cpp`: snapshot the iteration-mutable state
  (`isol`/`vsol` on the heap — `stagpoint_move`→`rebuild_ue_m` rewrites
  `vsol.ue_m`; `update_transition` rewrites `turb`/ctau — plus `glob.U` and the
  freeze history), try the limiter-ω step, accept on a residual criterion, else
  halve and roll back (≤5 backtracks). Gated on stalling (`r0 > 0.9·prev`) so
  clean cases take the unchanged damped-Newton step. All control logic plain
  `double` via `.getValue()`. Regression stayed 10/10 bit-identical.

## Why it does not work here

1. **The correct merit is the *settled* residual.** A first attempt evaluated
   the trial residual right after the Newton state update (before
   `stagpoint_move`/`update_transition`). In this viscous–inviscid coupled
   solver the residual at that half-finished configuration is not comparable to
   `r0`: the trace showed the post-step residual *above* `r0` on every iteration
   of the clean NACA 0012 case even though the *settled* residual decreased
   ~2 %/iter — so the search throttled every good step to ω/32 and the golden
   case went from 9 iterations to never converging. Fixing the merit to the
   settled residual (re-run stag-move + transition per trial, with rollback)
   restored bit-identical clean-case behaviour.

2. **Convergence here is intrinsically non-monotone.** Even with the settled
   merit, requiring strict descent (`trial ≤ r0·(1−c·ω)`) regressed the grid to
   78.0 % and tripled the `transition_front_oscillation` count (74 → 196). The
   solver routinely takes productive steps that *raise* the global residual
   while the transition front and stagnation point reorganise, then fall; strict
   descent rejects exactly these steps and traps the iteration in a slow crawl
   that hits the 60-iter cap.

3. **A non-monotone (Grippo-style) criterion recovers Phase 1 but adds nothing.**
   Accepting any trial within `1.2×` the recent residual envelope (rejecting only
   genuine blow-ups) returns to 86.0 % — statistically identical to Phase 1
   (Δ = −1 case) — while adding the snapshot/restore machinery and ~60 % runtime.

4. **The failures are not divergence.** After Phase 1, `diverged` is only 24/1375
   (1.7 %). The dominant non-converged buckets are limit cycles and stalls
   (`transition_front_oscillation`, `no_convergence`) where the residual
   oscillates in a bounded band with no descent trend — a line search cannot
   break a limit cycle, it can only slow the approach to it. The real lever is
   the transition-front / ctau dynamics (Phase 3), not step-length control.

## Decision

Line search reverted. Phase 1's RMS criterion is retained. Effort redirected to
**Phase 3 (ctau equilibrium seeding on transition movement)**, which targets the
`transition_front_oscillation` root cause that actually dominates the residual
failure set.

Artifacts: `sweep_p2.jsonl` (strict), `sweep_p2b.jsonl` (non-monotone) and their
logs are kept as evidence.
