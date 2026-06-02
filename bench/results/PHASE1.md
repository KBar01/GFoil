# Phase 1 — RMS convergence criterion

## Change

`coupled.cpp` previously converged on a raw, un-normalised Euclidean residual
norm over `Rsize = 3*(Ncoords+Nwake) ≈ 1500` entries against `rtol = 1e-10` —
roughly 5–6 decades tighter than XFOIL's RMS criterion, forcing the solver deep
into the regime where transition-front limit cycles dominate.

- `euc_norm` → `resid_rms`: divides the sum of squares by the entry count before
  the sqrt, so the tolerance is a size-independent per-equation RMS residual.
- `Param::rtol` default `1e-10 → 1e-6` (now interpreted as the RMS tolerance).
- `rtol` plumbed from `input["rtol"]` through `runCode` (default 1e-6) — a clean
  forward-only knob; lets the harness sweep tolerance without rebuilding. AD path
  untouched (`runCode`/`solve_coupled` are forward-TU only).

## Tolerance sweep (full 1375-case grid, cold-start)

| tolerance | cold-converged | median iters |
|---|---|---|
| baseline (1e-10 Euclidean) | 1175/1375 = 85.5 % | 16 |
| RMS 1e-5 | 1185 = 86.2 % | 12 |
| **RMS 1e-6 (chosen)** | **1184 = 86.1 %** | **13** |
| RMS 1e-7 | 1181 = 85.9 % | 13 |

Tolerances 1e-5…1e-7 are essentially flat on convergence. **1e-6** is chosen: it
keeps the NACA 0012 golden case bit-identical to the old 1e-10 result (converges
at the same iteration; CL/CD/CM/OASPL match to ~1e-11) while giving the iteration
speedup, whereas 1e-5 drifts the golden CL at the 7th digit.

## Effect

- **Speed**: median converged Newton iterations 16 → 13 (~20 % fewer); the golden
  case drops 10 → 9 iters, anchors drop 1–3 iters each.
- **Robustness**: +0.6 % cold-converged (+9 cases). Modest — as expected. The
  `diverged` bucket collapses 108 → 24, but most of those former-"diverged" cases
  are **reclassified**, not converged: they were never diverging (residual ≥ 1.0
  Euclidean ≈ RMS 0.026), they were stalling in a limit cycle at moderate
  residual. Under the RMS metric they now correctly read as
  `transition_front_oscillation` / `no_convergence`. This is more accurate
  labelling, and it isolates the real problem for Phase 2: genuine limit cycles
  and ctau-limiter stalls that oscillate at RMS ≫ 1e-6 and a looser tolerance
  cannot reach. The line search (Phase 2) is the lever for those.

## Anchors (rtol=1e-6) — no regressions

| anchor | cold-start | iters (baseline → P1) |
|---|---|---|
| NACA 0012 nCrit=5 α=±2.5° | CONVERGES | 7 → 6 |
| NACA 0012 nCrit=5 α=+4.7° | CONVERGES | 13 → 11 |
| NACA 0008-34 α=−2.6° | FAILS (osc) | 100 → 100 |
| B737 Midspan α=−3.1/−3.2° | CONVERGES | 36/35 → 34/33 |

## Golden

Regenerated. Forward scalars bit-identical to the previous golden (≤1e-11 rel);
AD gradients shifted ≤1.6e-5 rel because the converged linearisation point moves
by one Newton iteration (iter 9 vs 10). Not a correctness change — XFOIL
agreement on CL/CD/CM is preserved.
