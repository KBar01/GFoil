# Part B — Levenberg–Marquardt nan_lock rescue: investigated and rejected

## Result

LM diagonal regularisation of the singular forward Newton solve **rescued 0 of
the 24 `nan_lock` cases on the benchmark grid, with 0 regressions** — net zero.
Reverted. `src/` after Part B equals the post-A state.

| | post-A | post-B (LM) |
|---|---|---|
| full cold-converged | 1203/1375 = 87.5 % | 1203/1375 = 87.5 % |
| nan_lock | 24 | 24 |
| regressions | — | 0 |

## What was built

`solve_sys_sparse` (sparselinsolve.hpp): when the pure-Newton factorisation fails
or yields a non-finite step, retry `(A + λ·D) dU = b` with `D` = scale-aware
absolute diagonal of `A`, `λ` ramped 1e-6 → ×10 up to 6 times, falling back to
the existing `dU = 0` only if all attempts fail.

Adjoint safety was clean by construction (and confirmed): the forward solve runs
with an **inactive tape** (no `setActive` anywhere in the forward path), so the
external-function block is never built during the forward solve, and the AD
gradients come from the entirely separate `solve_sys_ad` path. The LM branch is
entered **only** when the pure-Newton solve fails, so the golden and every
currently-solvable case are bit-identical. Regression stayed **10/10**.

## Why it does not work here

`GFOIL_DEBUG` traces on the benchmark (smoothed) geometries show
`[LM] regularisation failed (A_finite=0)` for every `nan_lock` case: **the BL
Jacobian contains NaN/Inf entries, not finite-singular ones.** Diagonal
regularisation is structurally powerless against this — `NaN + λ·D = NaN`, so
`A + λD` never factorises. This is not an ill-conditioning problem that a
trust-region damps; it is a genuine rank-deficient / NaN-producing boundary-layer
state (the `get_*` closures returning non-finite derivatives at a near-critical
node). As the task anticipated, "a singular Jacobian that stays singular under
diagonal regularisation indicates a genuine rank-deficient BL state (true
attractor) … warm-start continuation remains the answer."

(A transient apparent success of "4/24" during development was a test-harness
geometry bug — the focused probe loaded root-directory `.dat` files for foils
that also have a smoothed-folder copy, AG11 / NACA 0008-34 / USA-35B, i.e. a
different geometry than the benchmark. On the benchmark geometries the gain is 0,
as the full sweep confirms.)

## Decision

Reverted. The 24 `nan_lock` cases are NaN-valued-Jacobian attractors outside the
reach of linear-solve regularisation; they remain a warm-start-continuation
class. Part A's acoustic fix stands as the deliverable of this task.

Artifact: `sweep_postB.jsonl` + `sweep_postB.log`.
