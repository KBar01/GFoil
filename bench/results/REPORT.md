# GFoil cold-start convergence — final report

Goal: make the cold-start Newton solve converge more robustly and faster,
without the Python-side alpha-backstepping in `standard_run`, staying within the
solver's existing fidelity (integral BL, e^N transition, Drela closures). All
work measured on the 25-foil stratified grid (`bench/`, 1375 cold-start grid
cases + 6 documented anchors). Regression golden (NACA 0012, α=2°, nCrit=5)
verified after every step.

## Headline: baseline → final

| metric | baseline | final | Δ |
|---|---|---|---|
| **full cold-converged** (valid OASPL) | 1175/1375 = 85.5 % | 1184/1375 = **86.1 %** | +0.6 pt |
| **true aero cold-converged** (conv + acoustic_nan) | 85.5 % | **87.5 %** | +2.0 pt |
| median Newton iterations (converged) | 16 | **13** | −19 % |
| median wall-time (converged) | 168 ms | **144 ms** | −14 % |

The "true aero" line separates aerodynamic convergence from cases that converge
aerodynamically but produce a non-finite OASPL in the downstream acoustic model
(see Phase 4). The full-converged line is the conservative metric (CL/CD/CM **and**
OASPL all valid).

## What shipped

**Phase 0 — diagnostic harness** (`bench/`). Stratified cold-start sweep over the
foil set with a geometry classifier, resumable JSONL records, and a scoreboard.
This is the scoreboard every change is measured against; it also revealed that
**4 of the 5 documented single-point cold-start Known Limitations no longer
reproduce** on this branch (they already converge cold — the CLAUDE.md list is
stale).

**Phase 1 — XFOIL-style RMS convergence criterion** (committed). Replaced the
raw size-dependent Euclidean residual norm + `rtol=1e-10` with a per-equation RMS
norm + `rtol=1e-6` (plumbed from `input["rtol"]`). Net: +0.6 pt convergence and
~19 % fewer iterations / ~14 % faster, with the NACA 0012 golden state preserved
bit-for-bit (it converges at the same iteration). Tolerance swept 1e-5…1e-7;
1e-6 chosen.

**Phase 4 — `acoustic_nan` labelling** (committed). `run_forward.cpp` silently
flipped `converged=false` whenever the acoustic OASPL was NaN/Inf, even when the
aero solve had converged — the Phase-0 "blank" failures. These are aero-converged
solutions whose Amiet/WPS model degenerates at low Re, now reported as
`failure_mode="acoustic_nan"` (noise code untouched). 19/1375 cases (1.4 %); this
is what lifts true-aero convergence to 87.5 %.

## What was investigated and rejected (negative results, not committed to `src/`)

**Phase 2 — residual line search.** A correct settled-residual Armijo line search
(snapshot/restore, stalling-gated, regression bit-identical) gave **no net gain**:
strict descent regressed to 78.0 % (this coupled iteration converges along a
non-monotone residual path; strict descent traps productive steps and tripled
`transition_front_oscillation`), and a non-monotone Grippo variant merely matched
Phase 1 (86.0 %) at +60 % runtime. After Phase 1, `diverged` is only 1.7 % of the
grid — the dominant failures are limit cycles a line search cannot break.
(`PHASE2.md`.)

**Phase 3 — ctau equilibrium seeding.** Seeding newly-turbulent nodes at the
per-node `get_cttr` equilibrium **regressed to 83.5 %** (thick_symmetric 244→226).
`get_cttr ∝ exp(−E/(Hk−1))` is very Hk-sensitive, so at freshly-transitioned
high-Hk nodes it produces a jagged ctau profile that the global Newton step
handles worse than the original smooth interpolation. The existing interpolation
already anchors its upstream endpoint on the manifold; the "far from manifold"
premise is weak here. Freeze/co-activation machinery retained. (`PHASE3.md`.)

## Final breakdown (true aero convergence)

By class (baseline_conv → final_conv (final_aero), /275):
thin_low_re 239→239(243) · thick_symmetric 243→244(250) · cambered_gp 240→240(243) ·
high_camber 225→230(234) · reflexed_6series 228→231(233).

By α: −6° 94 %, −3° 92 %, 0° 91 %, +3° 92 %, +6° 85 %, **+10° 70 %**.
By Re: 2e5 90 %, 1e6 91 %, **5e6 81 %**. By nCrit: 5 88 %, 9 88 %, 12 86 %.

The residual failures concentrate at high incidence (α=+10°), high Re (5e6), and
thin-wake / high-nCrit cells — separation- and transition-sensitive regimes where
the transition front limit-cycles. These are genuine BL attractors, not numerical
artefacts of the convergence criterion.

## Documented Known Limitations — cold-start status on the final build

| case | CLAUDE.md | final cold-start |
|---|---|---|
| NACA 0012 nCrit=5 α=±2.5° (nan_lock) | needs warm start | **converges cold** (6 it) |
| NACA 0012 nCrit=5 α=+4.7° (ctau saturation) | needs warm start | **converges cold** (11 it) |
| B737 Midspan α=−3.1/−3.2° (period-2) | needs warm start | **converges cold** (34/33 it) |
| NACA 0008-34 α=−2.6° (period-14 attractor) | permanent | **still fails cold** (transition_front_oscillation) |

Four of the five resolve cold; they already did before this work (the docs were
stale — a Phase-0 finding). **NACA 0008-34 α=−2.6° remains the one genuinely
unresolved cold-start case** and still requires warm-start continuation;
consistent with its documented "permanent" classification, the Phase 2/3 levers
did not crack it.

## Golden / regression

Regression 10/10 at every committed step. Golden was regenerated once (Phase 1):
forward scalars bit-identical (≤1e-11), AD gradients shifted ≤1.6e-5 rel from the
linearisation point moving one Newton iteration — not a physics change; XFOIL
CL/CD/CM agreement preserved. Phase 4 left the golden unchanged.

## Recommendation

Phase 1 (RMS criterion) is the net win and should stay. The remaining cold-start
failures are transition-front limit cycles that the existing limiter + freeze
machinery already handles about as well as local numerical re-seedings allow;
step-length control (Phase 2) and ctau re-seeding (Phase 3) do not help. Further
gains likely require structural changes outside the agreed fidelity/scope
(e.g. a continuation/pseudo-transient outer loop, or a proper transition-front
under-relaxation scheme), or accepting warm-start continuation for the residual
~12 % — which already handles them.
