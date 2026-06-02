# Phase 3 — ctau equilibrium seeding: investigated and rejected

## Summary

Seeding newly-turbulent nodes' ctau at the local `get_cttr` equilibrium (instead
of the existing linear interpolation) **regressed** cold-start convergence. It
was **reverted**; the `src/` tree after Phase 3 equals Phase 1.

| variant | cold-converged | median iters |
|---|---|---|
| Phase 1 | 1184/1375 = 86.1 % | 13 |
| Phase 3 per-node `get_cttr` seed | 1148/1375 = 83.5 % | 12 |

Per geometry class (P1 → P3): thick_symmetric **244 → 226**, high_camber
230 → 220, reflexed_6series 231 → 226, cambered_gp 240 → 237, thin_low_re 239 →
239. `no_convergence` rose 50 → 79.

## What was changed

`update_transition.cpp` (transition-advance branch): replaced

```cpp
sa0 = get_cttr(front node);  sa1 = current ctau at old first-turbulent node;
for each new node:  ctau = lerp(sa0, sa1, arc-fraction);   // smooth, off-manifold
```

with a per-node equilibrium seed

```cpp
for each new node:  ctau = get_cttr(that node's own th,ds,ue);  // on-manifold, jagged
```

Regression moved only at the 7th–8th significant figure (CL/CD/CM ≤1e-7,
gradients ≤5e-5) — the same physical fixed point via a different path — so the
change was correct and CoDi-safe; it simply converges *less* often.

## Why it does not work here

`get_cttr = CtauC·exp(−CtauE/(Hk−1))·cteq` is extremely sensitive to `Hk`. At a
just-transitioned node the boundary-layer shape (`th`, `ds`) still carries the
high-`Hk` laminar/transitional signature, so the per-node equilibrium ctau
varies sharply from node to node, producing a **jagged** ctau seed across the
newly-turbulent block. The global Newton step handles the original **smooth**
interpolated profile (front-equilibrium → existing turbulent ctau) better than
this jagged on-manifold profile, even though the latter is "closer to
equilibrium" pointwise. The interpolation is already a well-tuned compromise:
its upstream endpoint is the `get_cttr` equilibrium and its downstream endpoint
is the (near-converged, hence near-equilibrium) existing turbulent ctau, with a
smooth ramp between — so the premise that the interpolated seed is "far from the
manifold" does not hold strongly in practice.

A refinement anchoring both interpolation endpoints on the manifold (replace the
downstream endpoint with its `get_cttr`) was considered but rejected without a
sweep: the old first-turbulent node is already near-converged, so its current
ctau ≈ its equilibrium, making that variant ≈ the original.

## Interaction with the freeze machinery

The task suggested simplifying/removing the `ctau_freeze` / co-activation logic
*once seeding is correct*. Since improved seeding did not materialise (it
regressed), removing the freeze machinery — which demonstrably suppresses the
period-N transition-front cycles in `coupled.cpp` — would only make things
worse. The freeze logic is retained unchanged.

## Decision

Reverted. Combined with Phase 2, the conclusion is that this solver's cold-start
failures are dominated by transition-front limit cycles that the existing
limiter + freeze machinery already handles about as well as these local
numerical re-seedings allow. Phase 1's RMS criterion remains the net
improvement. See the final report (`bench/results/REPORT.md`).

Artifact: `sweep_p3.jsonl` + `sweep_p3.log`.
