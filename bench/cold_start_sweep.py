#!/usr/bin/env python3
"""Cold-start convergence sweep over the Smoothed_TEfixed_linear geometry library.

Samples 60 geometries (stratified over the sorted filename list: the list is
split into 60 equal contiguous strata and one file is drawn per stratum with
numpy default_rng(seed=0)), then runs the cold-start grid
alpha in [-4,-2,0,2,4,6,8] x nCrit in {5,9} at Re=2e6, rtol=1e-6 — 840 solves.

"Cold start" = the single first solve that fwd_run/standard_run performs
(gfoil_cpp.run_forward with no restart state). For every cold failure the full
fwd_run continuation (alpha backstepping) is also run to record whether the
driver rescues the case (rescue_conv column).

Each foil runs in its own subprocess so a hard crash (e.g. malformed geometry)
costs only that foil's remaining cases (failure_mode="crashed").

Rerun:
  source /home/pa20830/NewGradientVal/venv/bin/activate
  python3 bench/cold_start_sweep.py            # sweep + report
  python3 bench/cold_start_sweep.py --report   # rebuild REPORT.md from the CSV

Conditions held at defaults: Ma=0, rho=1.225, nu=1.5e-5, chord=1, span=3
(Aerofoil default), panelUniformity=1.0, panelTEspacing=0.09, free transition,
model "kam", TESampleLoc=0.98 (scalar), observer (0,3,0.5), aWeighting off,
repanel=False.
"""

import argparse
import csv
import json
import math
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FOIL_DIR = REPO / "Smoothed_TEfixed_linear"
RESULTS_DIR = REPO / "bench" / "results"
CSV_PATH = RESULTS_DIR / "cold_start_sweep.csv"
REPORT_PATH = RESULTS_DIR / "REPORT.md"

ALPHAS = [-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0]
NCRITS = [5.0, 9.0]
N_FOILS = 60
SEED = 0
RE = 2e6
RTOL = 1e-6

CSV_FIELDS = ["foil", "alpha", "ncrit", "conv", "aero_conv", "failure_mode",
              "newton_iterations", "wall_ms", "CL", "CD", "OASPL", "rescue_conv"]


def sample_foils():
    files = sorted(p.name for p in FOIL_DIR.glob("*.dat"))
    if len(files) < N_FOILS:
        raise SystemExit(f"only {len(files)} .dat files in {FOIL_DIR}")
    import numpy as np
    rng = np.random.default_rng(SEED)
    picks = []
    for s in range(N_FOILS):
        lo = s * len(files) // N_FOILS
        hi = (s + 1) * len(files) // N_FOILS
        picks.append(files[int(rng.integers(lo, hi))])
    return picks


def load_dat(path: Path):
    """Selig/Tecplot-style .dat: skip non-numeric lines, return (x, y) lists."""
    xs, ys = [], []
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        xs.append(x)
        ys.append(y)
    return xs, ys


# ---------------------------------------------------------------------------
# Worker: one foil, all grid cases, JSONL out
# ---------------------------------------------------------------------------

def run_worker(foil_name: str, out_path: Path) -> None:
    sys.path.insert(0, str(REPO))
    import numpy as np
    from GFoil.inputs import Aerofoil, OperatingConds, Acoustics
    from GFoil.gfoil import _build_input_dict, fwd_run
    import GFoil.gfoil_cpp as cpp

    xs, ys = load_dat(FOIL_DIR / foil_name)
    foil = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0, 3.0, 0.5]]))  # model/TESampleLoc defaults

    with open(out_path, "w") as fh:
        for ncrit in NCRITS:
            for alpha in ALPHAS:
                op = OperatingConds(alpha=alpha, Re=RE, nCrit=ncrit, rtol=RTOL)
                inp = _build_input_dict(foil, op, ac)
                t0 = time.perf_counter()
                r = cpp.run_forward(inp)        # the cold attempt
                wall_ms = (time.perf_counter() - t0) * 1e3
                conv = int(r.get("conv", 0))
                fm = r.get("failure_mode", "") or ""
                row = {
                    "foil": foil_name, "alpha": alpha, "ncrit": ncrit,
                    "conv": conv,
                    "aero_conv": int(conv == 1 or fm == "acoustic_nan"),
                    "failure_mode": fm if not conv else "",
                    "newton_iterations": r.get("newton_iterations", -1),
                    "wall_ms": round(wall_ms, 1),
                    "CL": r.get("CL") if conv else None,
                    "CD": r.get("CD") if conv else None,
                    "OASPL": r.get("OASPL") if conv else None,
                    "rescue_conv": "",
                }
                if not conv:
                    rr = fwd_run(foil, op, ac)   # full driver with backstepping
                    row["rescue_conv"] = int(rr.converged)
                fh.write(json.dumps(row) + "\n")
                fh.flush()


# ---------------------------------------------------------------------------
# Driver: subprocess per foil
# ---------------------------------------------------------------------------

def run_sweep() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    foils = sample_foils()
    rows = []
    t_start = time.perf_counter()
    for i, foil in enumerate(foils, 1):
        tmp = RESULTS_DIR / ".worker.jsonl"
        if tmp.exists():
            tmp.unlink()
        p = subprocess.run(
            [sys.executable, __file__, "--worker", foil, "--worker-out", str(tmp)],
            capture_output=True, text=True)
        done = []
        if tmp.exists():
            done = [json.loads(l) for l in tmp.read_text().splitlines() if l.strip()]
        rows += done
        seen = {(d["alpha"], d["ncrit"]) for d in done}
        for ncrit in NCRITS:                     # fill rows lost to a crash
            for alpha in ALPHAS:
                if (alpha, ncrit) not in seen:
                    rows.append({"foil": foil, "alpha": alpha, "ncrit": ncrit,
                                 "conv": 0, "aero_conv": 0,
                                 "failure_mode": "crashed",
                                 "newton_iterations": -1, "wall_ms": -1,
                                 "CL": None, "CD": None, "OASPL": None,
                                 "rescue_conv": ""})
        n_ok = sum(d["conv"] for d in done)
        crash = "" if p.returncode == 0 else f"  [worker exit {p.returncode}]"
        print(f"[{i:2d}/{len(foils)}] {foil:<45s} {n_ok:2d}/14 cold-converged{crash}")
    if (RESULTS_DIR / ".worker.jsonl").exists():
        (RESULTS_DIR / ".worker.jsonl").unlink()

    with open(CSV_PATH, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} cases -> {CSV_PATH} "
          f"({(time.perf_counter()-t_start)/60:.1f} min)")


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def pct(a, b):
    return f"{100.0*a/b:.1f}%" if b else "n/a"


def make_report() -> None:
    with open(CSV_PATH) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["alpha"] = float(r["alpha"]); r["ncrit"] = float(r["ncrit"])
        r["conv"] = int(r["conv"]); r["aero_conv"] = int(r["aero_conv"])
        r["newton_iterations"] = int(r["newton_iterations"])
        r["wall_ms"] = float(r["wall_ms"])

    n = len(rows)
    conv = [r for r in rows if r["conv"]]
    aero = [r for r in rows if r["aero_conv"]]
    fails = [r for r in rows if not r["conv"]]

    L = []
    L.append("# Cold-start convergence baseline — 60-foil sweep\n")
    L.append(f"Library: `Smoothed_TEfixed_linear/` ({len(sorted(FOIL_DIR.glob('*.dat')))} files); "
             f"60 foils stratified-sampled (seed={SEED}, sorted order, one per contiguous stratum).")
    L.append(f"Grid: alpha {ALPHAS} x nCrit {[int(x) for x in NCRITS]}, Re={RE:g}, rtol={RTOL:g}, "
             "repanel=False, free transition, model kam, TESampleLoc 0.98, observer (0,3,0.5), span 3.")
    L.append("Cold start = single `gfoil_cpp.run_forward` call (what fwd_run tries first); "
             "`rescue_conv` = full fwd_run backstepping outcome for cold failures.\n")
    L.append("Rerun: `python3 bench/cold_start_sweep.py` "
             "(venv: /home/pa20830/NewGradientVal/venv)\n")

    L.append("## Headline\n")
    L.append(f"- cold-converged (valid OASPL): **{sum(r['conv'] for r in rows)}/{n} = {pct(len(conv), n)}**")
    L.append(f"- true aero cold-converged (conv + acoustic_nan): **{len(aero)}/{n} = {pct(len(aero), n)}**")
    resc = [r for r in fails if r["rescue_conv"] == "1"]
    L.append(f"- cold failures rescued by fwd_run backstepping: {len(resc)}/{len(fails)}")
    its = sorted(r["newton_iterations"] for r in conv)
    if its:
        p90 = its[int(0.9 * (len(its) - 1))]
        L.append(f"- Newton iterations (converged): median {int(statistics.median(its))}, "
                 f"p90 {p90}, max {its[-1]}")
        ms = sorted(r["wall_ms"] for r in conv)
        L.append(f"- wall ms (converged): median {statistics.median(ms):.0f}, "
                 f"p90 {ms[int(0.9*(len(ms)-1))]:.0f}, max {ms[-1]:.0f}\n")

    L.append("## Convergence by alpha\n")
    L.append("| alpha | cold conv | true aero |")
    L.append("|---|---|---|")
    for a in ALPHAS:
        g = [r for r in rows if r["alpha"] == a]
        L.append(f"| {a:+.0f} | {pct(sum(r['conv'] for r in g), len(g))} "
                 f"| {pct(sum(r['aero_conv'] for r in g), len(g))} |")

    L.append("\n## Convergence by nCrit\n")
    L.append("| nCrit | cold conv | true aero |")
    L.append("|---|---|---|")
    for nc in NCRITS:
        g = [r for r in rows if r["ncrit"] == nc]
        L.append(f"| {nc:g} | {pct(sum(r['conv'] for r in g), len(g))} "
                 f"| {pct(sum(r['aero_conv'] for r in g), len(g))} |")

    L.append("\n## Failure modes\n")
    L.append("| failure_mode | count | % of grid | rescued by backstepping |")
    L.append("|---|---|---|---|")
    modes = {}
    for r in fails:
        modes.setdefault(r["failure_mode"], []).append(r)
    for m, g in sorted(modes.items(), key=lambda kv: -len(kv[1])):
        nresc = sum(1 for r in g if r["rescue_conv"] == "1")
        L.append(f"| {m} | {len(g)} | {pct(len(g), n)} | {nresc}/{len(g)} |")

    L.append("\n## Slowest-converging 10 cases\n")
    L.append("| foil | alpha | nCrit | iters | ms |")
    L.append("|---|---|---|---|---|")
    for r in sorted(conv, key=lambda r: -r["newton_iterations"])[:10]:
        L.append(f"| {r['foil']} | {r['alpha']:+.0f} | {r['ncrit']:g} "
                 f"| {r['newton_iterations']} | {r['wall_ms']:.0f} |")

    L.append("\n## All failing cases (file, alpha, nCrit, mode, rescued)\n")
    for r in sorted(fails, key=lambda r: (r["failure_mode"], r["foil"], r["alpha"])):
        L.append(f"- `{r['foil']}` alpha={r['alpha']:+.0f} nCrit={r['ncrit']:g} "
                 f"— {r['failure_mode']}"
                 + (" (rescued)" if r["rescue_conv"] == "1" else ""))

    REPORT_PATH.write_text("\n".join(L) + "\n")
    print(f"report -> {REPORT_PATH}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", default=None, help=argparse.SUPPRESS)
    ap.add_argument("--worker-out", default=None, help=argparse.SUPPRESS)
    ap.add_argument("--report", action="store_true",
                    help="only rebuild REPORT.md from the existing CSV")
    args = ap.parse_args()

    if args.worker:
        run_worker(args.worker, Path(args.worker_out))
    elif args.report:
        make_report()
    else:
        run_sweep()
        make_report()


if __name__ == "__main__":
    main()
