# GFoil nomenclature

Single reference for the terse, domain-standard symbol names used throughout the
solver — especially in the boundary-layer closure relations (`get_funcs.hpp`),
residuals (`residuals_shared.hpp`), and the Newton solve. These names are the
**lingua franca of integral boundary-layer methods** (Drela's XFOIL/MSES, and the
Fidkowski (2021) formulation this code follows); they are deliberately **kept**
rather than renamed so the code reads against the literature. This glossary
defines each one once.

References cited in-code: *Fidkowski, "A Coupled Inviscid–Viscous Airfoil Analysis
Solver, Revisited" (2021)* (equation numbers in comments), and Drela & Giles,
*"Viscous–Inviscid Analysis of Transonic and Low Reynolds Number Airfoils"* (AIAA
J. 1987), the origin of the e^N transition and turbulent closure relations.

## Primary BL state vector `U`

Each surface node carries a 4-component state. Throughout the code the per-node
state is indexed `U[colMajorIndex(k, node, 4)]` with `k`:

| k | symbol | meaning | units |
|---|--------|---------|-------|
| 0 | `th`   | momentum thickness θ | length |
| 1 | `ds`   | displacement thickness δ* | length |
| 2 | `sa`   | **dual-use**: laminar → amplification factor ñ (e^N); turbulent → shear-stress coefficient c_τ | – |
| 3 | `ue`   | boundary-layer edge velocity (signed; sign encodes surface side) | velocity |

The `_U` suffix on a derivative output (e.g. `Hk_U`, `cf_U`) means **the gradient
of that quantity with respect to this 4-vector** `[∂/∂th, ∂/∂ds, ∂/∂sa, ∂/∂ue]`.
A trailing `_x` (e.g. `cfxt_x`) is the derivative w.r.t. the streamwise arc-length
coordinate ξ.

## Shape factors and BL quantities

| symbol | meaning |
|--------|---------|
| `H`    | shape factor δ*/θ (`get_H`) |
| `Hk`   | kinematic shape factor — `H` with the Karman–Tsien compressibility correction removed (`get_Hk`) |
| `Hs`   | kinetic-energy shape factor H* = θ*/θ (`get_Hs`) |
| `Hss`  | density shape factor H** (`get_Hss`) |
| `Hw`   | wake-gap shape factor (wake gap / θ) (`get_Hw`) |
| `de`   | boundary-layer thickness δ (the 99% thickness; `get_de`), capped at 12·θ |
| `Ret`  | momentum-thickness Reynolds number Re_θ = ρ·θ·u_e/μ (`get_Ret`) |
| `Us`   | normalised slip velocity (Drela's U_s) (`get_Us`) |
| `uq`   | equilibrium-locus parameter U_q feeding the lag (`get_uq`) |

## Friction, dissipation, shear

| symbol | meaning |
|--------|---------|
| `cf`   | skin-friction coefficient C_f (`get_cf`; laminar Fidkowski Eq.56, turbulent Eqs.57–60) |
| `cfxt` | C_f scaled into the momentum-equation residual term, C_f·ξ/θ (`get_cfxt`) |
| `cDi`  | dissipation coefficient C_D (`get_cDi` and the `_lam/_lamstress/_lamwake/_outer/_turbwall` component builders) |
| `cDixt`| C_D scaled into the shape-parameter residual term (`get_cDixt`) |
| `cteq` | equilibrium shear-stress coefficient c_{τ,eq} (`get_cteq`) |
| `cttr` | shear-stress coefficient seeded at transition (`get_cttr`) |
| `tauWall` / `tauMax` | wall / maximum shear stress [Pa] (acoustic TE inputs) |

## Compressibility / thermodynamics

| symbol | meaning |
|--------|---------|
| `uk`   | compressible edge speed from the incompressible one (Karman–Tsien, `get_uk`) |
| `Mach2` / `M2` | local Mach number squared (`get_Mach2`) |
| `cp`   | pressure coefficient C_p (`get_cp`) |
| `Minf`, `Vinf`, `Ma` | freestream Mach, freestream speed, Mach number |
| `rho0`, `mu0`, `H0`  | stagnation density, stagnation viscosity, stagnation total enthalpy |
| `Fc`   | compressibility factor on C_f (√(1+½(γ−1)M²)) |

## Transition / amplification

| symbol | meaning |
|--------|---------|
| `sa` (laminar) | envelope amplification factor ñ; transition where ñ reaches `ncrit` |
| `damp` | amplification rate dñ/dξ along the surface (`get_damp`; e^N envelope, Drela) |
| `ncrit` | critical amplification factor (transition threshold) |
| `ilam` | index of the last laminar node on a surface (transition front) |

## Discretisation / solve

| symbol | meaning |
|--------|---------|
| `upw`  | upwinding fraction for the two-point BL difference scheme (`get_upw`) |
| `xi` / `dist` / `distFromStag` | streamwise arc-length coordinate ξ measured from the stagnation point |
| `Is[si]` | ordered node-index list for surface `si` (0 = lower, 1 = upper, 2 = wake), stag → TE |
| `glob.R`, `glob.R_V_*` | global residual vector and its sparse Jacobian triplets |
| `glob.U`, `glob.dU` | global state and Newton update |
| `omega` | Newton relaxation / under-relaxation factor |
| `rtol` | RMS residual convergence tolerance (XFOIL-style; see `coupled.cpp`) |

## A note on the `Real` type

`Real` is `codi::RealReverse` (reverse-mode AD) in the active build and `double`
in passive contexts. Any arithmetic on `Real` inside an active region is recorded
on the CoDiPack tape, so **re-expressing or re-ordering it changes the adjoint**.
See the CoDiPack rules in `CLAUDE.md`. The `.getValue()` accessor extracts the
passive `double` and is only used for control-flow decisions and output.
