# GFoil cold-start convergence — Phase 0 baseline

Scoreboard against which every subsequent solver change is measured. Generated
by `bench/cold_start_sweep.py` on the `forcing-transition` branch, commit prior
to any Phase 1+ solver change.

## Method

- **Cold-start path only**: each case calls the in-process pybind11 solver via
  `GFoil.gfoil._call_forward(inp)` with `restart=0` and `prev_result=None`, i.e.
  bypassing `standard_run`'s alpha-backstepping/warm-start continuation.
- **Subset**: 25 foils, stratified 5-per-class over 5 geometry classes
  (`bench/results/subset.json`), selected automatically by parsed name + computed
  max thickness / camber / TE-reflex.
- **Grid per foil**: alpha ∈ {−6,−3,0,3,6,10}°, Re ∈ {2e5,1e6,5e6},
  nCrit ∈ {5,9,12} (54 cases) + 1 forced-transition case (x/c=0.3 both surfaces,
  alpha=3°, Re=1e6, nCrit=9). ncrithyst=0 throughout. Single observer, model=roz.
- **Anchors**: 6 fixed cases for the documented Known Limitations.
- 1375 grid + 6 anchor = **1381 cases**, total wall-time **337 s** (in-process).

## Headline

| metric | baseline |
|---|---|
| **Grid cold-converged** | **1175 / 1375 = 85.5 %** |
| Converged Newton iters | mean 20.0, median 16, max 59 |
| Converged wall-time | mean 205 ms, median 168 ms |

### Failure-mode breakdown (200 non-converged grid cases)

| mode | count | % of grid |
|---|---|---|
| diverged | 108 | 7.9 % |
| transition_front_oscillation | 28 | 2.0 % |
| nan_lock | 25 | 1.8 % |
| no_convergence | 21 | 1.5 % |
| (blank) | 18 | 1.3 % |

`diverged` (residual ≥ 1.0 at exit) is by far the largest bucket — the primary
target for Phase 1 (RMS criterion) and Phase 2 (line search). The 18 blank-mode
failures are all at Re=2e5 with no Python exception and varied iteration counts;
they appear to be early/pre-Newton exits where `failure_mode` is not propagated —
flagged for investigation but low volume.

### By geometry class

| class | converged |
|---|---|
| thin_low_re | 239/275 = 86.9 % |
| thick_symmetric | 243/275 = 88.4 % |
| cambered_gp | 240/275 = 87.3 % |
| high_camber | 225/275 = 81.8 % |
| reflexed_6series | 228/275 = 82.9 % |

### By operating cell (worst cells dominate the failures)

| alpha | conv | | Re | conv | | nCrit | conv |
|---|---|---|---|---|---|---|---|
| −6 | 91.1 % | | 2e5 | 86.2 % | | 5 | 86.0 % |
| −3 | 89.3 % | | 1e6 | 90.2 % | | 9 | 87.8 % |
| 0 | 87.1 % | | 5e6 | 79.3 % | | 12 | 82.0 % |
| +3 | 91.6 % | | | | | | |
| +6 | 84.0 % | | | | | | |
| **+10** | **68.4 %** | | | | | | |

High-incidence (α=+10°), high-Re (5e6), thin-wake / high-nCrit cells concentrate
the failures, consistent with separation- and transition-sensitivity.

Forced-transition cases: 24/25 = 96.0 % (forcing the front removes the
transition-tracking limit cycle, as expected).

## Documented anchors — current cold-start behaviour vs CLAUDE.md

| anchor | doc'd limitation | cold-start now |
|---|---|---|
| NACA 0012 nCrit=5 α=+2.5° | nan_lock | **CONVERGES** (7 it) |
| NACA 0012 nCrit=5 α=−2.5° | nan_lock | **CONVERGES** (7 it) |
| NACA 0012 nCrit=5 α=+4.7° | ctau saturation | **CONVERGES** (13 it) |
| NACA 0008-34 α=−2.6° | period-14 attractor (permanent) | FAILS (transition_front_oscillation) |
| B737 Midspan α=−3.1° | period-2 oscillation | **CONVERGES** (36 it) |
| B737 Midspan α=−3.2° | period-2 oscillation | **CONVERGES** (35 it) |

**Finding:** On this `forcing-transition` branch, 4 of the 5 documented
single-point cold-start limitations no longer reproduce — they already converge
cold. Only **NACA 0008-34 α=−2.6°** still fails cold, matching its
"permanent" classification. The CLAUDE.md Known Limitations list is stale
relative to current code and should be revisited once Phases 1–3 land.

## Files

- `bench/foil_select.py` — loader / geometry / classifier / subset selection
- `bench/cold_start_sweep.py` — the sweep + scoreboard (`--run`, `--report`)
- `bench/results/baseline.jsonl` — raw per-case records (resumable stream)
- `bench/results/subset.json` — the 25 selected foils with descriptors
- `bench/results/sweep_baseline.log` — run log (wall-time, per-failure trace)

Reproduce: `python3 bench/cold_start_sweep.py --run --fresh`
