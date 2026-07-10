# AD formulation audit — adjoint gradient correctness (July 2026)

Full audit of the mathematical formulation behind `grad_run`'s gradients
(`dCL_dy`, `dCD_dy`, `dOASPL_dy`, and the three alpha scalars): derivation,
code mapping, line-level verification of every hand-derived derivative, and
fresh numerical evidence. Run on branch `variableLength_input_geom` (post
variable-length-input change), 301-node `tests/input.json` conditions unless
stated.

**Verdict: the formulation is mathematically correct and the returned
gradients are exact derivatives of the forward solver as implemented.**
Strongest evidence: smooth-mode directional derivatives match central FD to
1.6e-7 … 6.9e-5 relative across free transition, forced transition, and
Ma=0.2, for all three outputs (table below). Two latent items to know about
(neither affects current gradients): `param.Minf` is never set so the
Karman–Tsien compressibility branches are dead code (finding F2), and
single-node FD checks near transition-front node flips are kinked at the
±1e-5 scale (finding F3 — an FD artefact, not an AD error).

---

## 1. The formulation the code implements

For output `g ∈ {CL, CD, OASPL}`, design variables `x = (y₁…y_n, α)`, and the
converged IBL state `U ∈ R^920` satisfying `R(U, x) = 0`:

```
dg/dx = ∂g/∂x − λᵀ ∂R/∂x,   with   (∂R/∂U)ᵀ λ = (∂g/∂U)ᵀ        (†)
```

Implemented in `run_AD_py` (`gfoil_ad_bindings.cpp`) in four steps:

1. **Pass 1 — `partialOutputspartialInputs<RealReverseVec<2>>`**
   (`ADfuncs.hpp`): tapes the *output evaluation only* (make_panels → foil →
   `calc_force` → `interpolate_at_95_both_surfaces` → `calc_OASPL`) with
   `y`, `α`, and all 920 states registered as independent inputs. Reverse
   sweep with the 2-vector seeded on (CL, OASPL) yields `∂g/∂y`, `∂g/∂α`,
   and `∂g/∂U` for both outputs in one evaluation.
2. **CD row — manual chain (`Realfwd`)**: CD is Squire–Young from the last
   wake node only, `CD = 2θ_w·ue_w^E`, `E = (5+H_w)/2`, `H_w = δ*_w/θ_w`
   (`calc_force`, solver_funcs.hpp; `Vinf = 1` in the nondim). `∂CD/∂U` has
   exactly three nonzero entries (θ, δ*, ue at node 229 = indices 916, 917,
   919), computed by hand with `codi::RealForward`. Verified symbolically
   here and numerically at Ma=0 (§3): the three entries match
   d(2θ·ue^E)/d(θ,δ*,ue) exactly. `∂CD/∂y = ∂CD/∂α = 0` is correct — the
   formula touches no geometry and no α explicitly.
3. **Adjoint solve — `solve_sys_ad`**: solves `(∂R/∂U)ᵀλ_g = ∂g/∂U` for the
   three outputs by SparseLU. The transpose is realised by swapping
   rows↔cols when unpacking the forward-saved triplets. Verified end to end:
   the forward stores `(row = residual index, col = state index)` in BOTH
   assembly paths (`equate_block_inplace_sparse`, build_global_sys.cpp:88,
   and the ue-coupling fill, solve_glob.cpp:66) and the unpack+call chain
   applies exactly one net transpose. Duplicate triplets are summed by
   `setFromTriplets` — correct semantics for the scatter-assembled Jacobian.
4. **Pass 2 — `partialRpartialx<RealReverseVec<3>>`**: tapes
   `h_g = −λ_gᵀ R(U, x)` with only `y` and `α` registered as inputs and `U`
   held passive (plain doubles), rebuilding the full residual (inviscid
   gammas → wake → stag → `build_glob_RV_AD` → `finishdRdU_AD`) from taped
   geometry. `∇h_g = −λ_gᵀ ∂R/∂x`, added to the pass-1 partials — exactly (†).

Alpha bookkeeping is consistent with (†): CL and OASPL get both terms
(`∂CL/∂α` from panel-force rotation; `∂OASPL/∂α` from the taped observer
rotation), CD gets the adjoint term only.

## 2. Line-level checks (all passed)

| Item | Verification |
|---|---|
| Adjoint identity, signs | `h = −Σλ_iR_i`; total = pass1 + pass2. Matches (†). |
| Jacobian transpose | Single net transpose confirmed by tracing the (deliberately confusing) double swap: load swap (`dRdU_rows[i]=RVcols[i]`, `dRdU_cols[i]=RVrows[i]`) into a `solve_sys_ad(cols, rows, …)` signature whose triplets are `(rows[k], cols[k])` → builds `A[state,residual] = (∂R/∂U)ᵀ`. Correct. |
| Jacobian evaluation point | `coupled.cpp` assembles `build_glob_RV` at the top of each iteration, tests convergence on that same residual, and saves that assembly (plus the ue-row Jacobian added by `solve_glob(...,0)`) — i.e. ∂R/∂U at the **accepted converged state**, all 920 rows. Consistent. |
| Residual consistency fwd↔AD | Both TUs call the *same* kernels: `residuals_shared.hpp` (`residual_station`, `residual_transition`, `wake_sys`) and `ue_residual_kernel` (solver_funcs.hpp). Forward instantiates `<true>` (analytic `R_U` blocks), AD instantiates `<false>` (residual only, CoDi-taped). So pass 2 differentiates the same `R` whose `∂R/∂U` the forward saved. |
| ue-coupling row Jacobian | `R_ue = ue − ue_inv − ue_m·(δ*∘ue)`; analytic fill `δ_rc − ue_m·δ*[c]` / `−ue_m·ue[c]` (solve_glob.cpp) — exact product rule. |
| Analytic BL `R_U` blocks | Not re-derived line-by-line here; covered by the June 2026 column-FD audit of dR/dU (CHANGELOG, DRDU audit: transition rows exact, stiff θ/δ* columns clean ∝h² truncation) and transitively by every gradient-level FD result below. |
| `errFunc` external fn | Pushed Jacobian is Cauchy–Riemann applied to `erf′(z) = 2/√π·e^{−z²}`: `u_x = v_y = Re w′`, `u_y = −Im w′`, `v_x = Im w′`. Correct. |
| `solve_sys_ue` external fn | For `AX = −B`: `Λ = A⁻ᵀX̄`, `B̄ −= Λ`, `Ā −= ΛXᵀ` — the textbook implicit-solve adjoint. The unregistered last row of X correctly carries zero adjoint. |
| `compute_Dw` external fn | For `Dw = Cgam·Bp + Csig`: `dCsig += L`, `dCgam += L·Bpᵀ`, `dBp += Cgamᵀ·L` — exact bilinear adjoint, linearised at recorded values (exact for a bilinear op). |
| Passive/active boundaries | `turb` pattern, stag bracket (`currStag`), transition-interval selection: integer/branch choices held passive (documented pattern); positions (`xift`, `distFromStag`, stag interpolation) taped in `Real`. `.getValue()` appears only in branch predicates and the sanctioned external-function/WriteJSON sites. |
| Tape hygiene | Separate tapes for Vec2/Vec3 types; inputs registered before use; outputs registered before `setPassive`; seed `gradient()[i]=1` per output; `reset()` after each pass + safety-net resets at entry. |

## 3. Numerical evidence (all runs `rtol = 1e-11`, central differences)

### Directional derivatives — smooth modes, whole-gradient-vector test

Modes: thickness `d_i = ±√x(1−x)` (sign per surface), camber `d_i = x(1−x)`.
`h = 1e-5` on `y + t·d`. Compares `∇g·d` (AD) vs FD:

| Case | dCL·d | dCD·d | dOASPL·d |
|---|---|---|---|
| free, Ma=0, thickness | 3.3e-05 | 1.2e-06 | 6.5e-07 |
| free, Ma=0, camber | 3.3e-07 | 2.6e-06 | 1.3e-06 |
| forced x/c=0.1, thickness | 1.5e-05 | 6.7e-08 | 6.1e-07 |
| forced x/c=0.1, camber | 1.6e-07 | 5.1e-07 | 1.2e-06 |
| free, Ma=0.2, thickness | 6.9e-05 | 1.3e-06 | 6.4e-07 |
| free, Ma=0.2, camber | 3.3e-07 | 3.0e-06 | 1.4e-06 |

(entries are relative errors)

### Alpha gradients (`h = 1e-4°`)

| | Ma=0 | Ma=0.2 |
|---|---|---|
| dCL/dα | 5.7e-08 | 6.0e-08 |
| dCD/dα | 1.3e-06 | 1.4e-06 |
| dOASPL/dα | 4.5e-07 | 5.0e-07 |

### Single-node h-sweeps (finding F3 illustration)

Node 30 (x = 0.905 lower, sharp-TE foil, Ma=0), AD `dCL/dy = 0.34084229`:
FD = 0.2020 (h=1e-4), 0.2941 (1e-5), 0.34080 (1e-6), 0.34089 (1e-7) —
**FD converges to the AD value**; the large-h values straddle a
transition-node-flip kink. Same signature reproduced on the 101-node coarse
case (node 35, one-sided slopes +0.11/−0.73) and at Ma=0.2. The regression
suite's FD spot check uses nodes screened away from the front (see
tests/README.md group 4).

## 4. Findings

**F1 — Formulation correct.** Every structural element of (†) is implemented
correctly; gradient values verified against FD across free/forced
transition, Ma=0/0.2, y-modes, y-nodes, and α.

**F2 — `param.Minf` is never assigned; Karman–Tsien is dead code (latent,
forward-side).** `Param::Minf` defaults to 0 and no line in src/ or srcAD/
sets it (`init_thermo` sets `KTb/KTl/H0/cps/mu0/rho0` but not `Minf`), so the
`param.Minf > 0` branches in `get_uk`, `get_Mach2`, `get_cp`, `get_Hk` can
never activate. At Ma=0.2 the forward's CD provably equals the raw-`ue`
formula (`2θ·ue^E = 0.00651349627` = reported CD; the KT value would be
0.00651111390). Ma currently affects the solve only through `mu0`/`rho0`.
Consequences: (a) the manual CD chain in the bindings — which differentiates
the raw-`ue` formula — is exactly consistent with the forward **as it
stands**, which is why all Ma=0.2 gradient checks pass; (b) if `Minf` is
ever wired up (`param.Minf = oper.Ma` in `init_thermo`), the manual chain
becomes inconsistent: at Ma=0.2 the ∂CD/∂δ* and ∂CD/∂ue entries would each
be ~2% wrong (`uk_ue = 1.0203` at the far-wake `ue = 0.9941`). Whoever
enables compressibility must update the `Realfwd` block in
`gfoil_ad_bindings.cpp` (and `srcAD/main.cpp`) to chain through
`uk(ue)` — and re-run the Ma>0 FD checks in §3.

**F3 — Transition-front kinks make naive per-node FD misleading (not a
bug).** The converged solution is C0-but-kinked in any design variable at
the points where the transition front crosses a panel node (the discrete
`turb` pattern flips; `residual_transition` moves to the adjacent interval).
AD returns the exact derivative of the locally-active smooth branch — the
correct object for optimisation, existing a.e. FD with h larger than the
distance to the nearest flip averages across branches and can read 10–100%
"error" (§3 h-sweeps) while FD(h→0) → AD everywhere tested. Verification
protocol: tight `rtol` (≤1e-10) AND either smooth modes or per-node h-sweeps.

**F4 — Fragility note (no action required).** The adjoint transpose is
spread across two swap sites with misleading names (`dRdU_rows` holds column
indices; `solve_sys_ad`'s first parameter is named `R_V_cols`). Net effect
is correct (§2), and `srcAD/main.cpp` replicates the same pattern
consistently, but any future edit touching either site should re-verify
against §3.

**F5 — Adjoint accuracy is bounded by forward convergence.** λ is computed
from ∂R/∂U and ∂g/∂U at the rtol-converged iterate, so gradient error scales
with the residual: at `rtol=1e-6` expect ~1e-4–1e-2 apparent FD error
(dominated by FD noise amplification 1/(2h), see CHANGELOG free-transition
FD artefact entry); at `rtol=1e-11` the agreement in §3 is recovered. This
is inherent to any adjoint on an inexactly-converged state, not a defect.

## 5. Reproduction

All experiments are ~20-line scripts against `GFoil.gfoil_cpp` following the
pattern of `tests/regression_test.py::run_coarse_fd_check` (base solve at
`rtol=1e-11` → `run_AD` → perturbed `run_forward` pairs). The directional
test perturbs all y simultaneously along the mode and compares `∇g·d`; the
permanent per-node spot check runs in the regression suite (48 checks).
