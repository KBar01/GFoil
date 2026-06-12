#!/usr/bin/env python3
"""Fast cold-only sweep — the 840 cold solves of the baseline grid, recording
per-case (conv, failure_mode, newton_iterations, CL). NO rescue/backstepping, so
this isolates the baseline contract (cold-converged count, per-alpha profile,
iteration distribution, failing-case set) and serves as the per-case trajectory
baseline for Part C.

Honors GFOIL_BFIX (passed through to each worker subprocess).

Usage:
  GFOIL_BFIX=<off|A|B> python3 bench/cold_only_sweep.py --out <csv>
  python3 bench/cold_only_sweep.py --compare <base.csv> <new.csv>
"""
import argparse
import csv
import json
import os
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
FOIL_DIR = REPO / "Smoothed_TEfixed_linear"

ALPHAS = [-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0]
NCRITS = [5.0, 9.0]
N_FOILS = 60
SEED = 0
RE = 2e6
RTOL = 1e-6
FIELDS = ["foil", "alpha", "ncrit", "conv", "failure_mode", "newton_iterations", "CL"]


def sample_foils():
    import numpy as np
    files = sorted(p.name for p in FOIL_DIR.glob("*.dat"))
    rng = np.random.default_rng(SEED)
    picks = []
    for s in range(N_FOILS):
        lo = s * len(files) // N_FOILS
        hi = (s + 1) * len(files) // N_FOILS
        picks.append(files[int(rng.integers(lo, hi))])
    return picks


def load_dat(path):
    xs, ys = [], []
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        xs.append(x); ys.append(y)
    return xs, ys


def run_worker(foil_name, out_path):
    import numpy as np
    from GFoil.inputs import Aerofoil, OperatingConds, Acoustics
    from GFoil.gfoil import _build_input_dict
    import GFoil.gfoil_cpp as cpp
    xs, ys = load_dat(FOIL_DIR / foil_name)
    foil = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0, 3.0, 0.5]]))
    with open(out_path, "w") as fh:
        for ncrit in NCRITS:
            for alpha in ALPHAS:
                op = OperatingConds(alpha=alpha, Re=RE, nCrit=ncrit, rtol=RTOL)
                inp = _build_input_dict(foil, op, ac)
                r = cpp.run_forward(inp)
                conv = int(r.get("conv", 0))
                fh.write(json.dumps({
                    "foil": foil_name, "alpha": alpha, "ncrit": ncrit,
                    "conv": conv,
                    "failure_mode": (r.get("failure_mode", "") or "") if not conv else "",
                    "newton_iterations": r.get("newton_iterations", -1),
                    "CL": r.get("CL") if conv else None,
                }) + "\n")
                fh.flush()


def run_sweep(out_csv):
    foils = sample_foils()
    rows = []
    for i, foil in enumerate(foils, 1):
        tmp = REPO / "bench" / "results" / ".cold_worker.jsonl"
        if tmp.exists():
            tmp.unlink()
        p = subprocess.run([sys.executable, __file__, "--worker", foil,
                            "--worker-out", str(tmp)],
                           capture_output=True, text=True, env=os.environ)
        done = []
        if tmp.exists():
            done = [json.loads(l) for l in tmp.read_text().splitlines() if l.strip()]
        rows += done
        n_ok = sum(d["conv"] for d in done)
        crash = "" if p.returncode == 0 else f" [exit {p.returncode}] {p.stderr[-200:]}"
        print(f"[{i:2d}/{len(foils)}] {foil:<45s} {n_ok:2d}/14{crash}")
        if tmp.exists():
            tmp.unlink()
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} cases -> {out_csv}")
    summarize(rows)


def summarize(rows):
    n = len(rows)
    conv = [r for r in rows if int(r["conv"])]
    print(f"cold-converged: {len(conv)}/{n}")
    print("by alpha:", {f"{float(a):+.0f}":
          sum(1 for r in rows if float(r["alpha"]) == a and int(r["conv"]))
          for a in ALPHAS})
    its = sorted(int(r["newton_iterations"]) for r in conv)
    if its:
        print(f"iters: median {int(statistics.median(its))}, "
              f"p90 {its[int(0.9*(len(its)-1))]}, max {its[-1]}")


def compare(base_csv, new_csv):
    base = {(r["foil"], r["alpha"], r["ncrit"]): r
            for r in csv.DictReader(open(base_csv))}
    new = {(r["foil"], r["alpha"], r["ncrit"]): r
           for r in csv.DictReader(open(new_csv))}
    keys = sorted(set(base) | set(new))
    conv_diff, fm_diff, it_diff = [], [], []
    over60 = []
    for k in keys:
        b, nw = base.get(k), new.get(k)
        if b is None or nw is None:
            print("MISSING", k); continue
        if b["conv"] != nw["conv"]:
            conv_diff.append((k, b["conv"], nw["conv"]))
        if (b.get("failure_mode", "") or "") != (nw.get("failure_mode", "") or ""):
            fm_diff.append((k, b.get("failure_mode"), nw.get("failure_mode")))
        if int(b["newton_iterations"]) != int(nw["newton_iterations"]):
            it_diff.append((k, int(b["newton_iterations"]), int(nw["newton_iterations"])))
        if int(nw["conv"]) and int(nw["newton_iterations"]) >= 60:
            over60.append((k, int(nw["newton_iterations"])))
    print(f"\n=== compare {Path(base_csv).name} -> {Path(new_csv).name} ===")
    print(f"conv-flag diffs:       {len(conv_diff)}")
    for k, a, b in conv_diff:
        print(f"   {k}  conv {a} -> {b}  (base fm={base[k].get('failure_mode')}, "
              f"new fm={new[k].get('failure_mode')})")
    print(f"failure_mode diffs:    {len(fm_diff)}")
    for k, a, b in fm_diff[:40]:
        print(f"   {k}  '{a}' -> '{b}'")
    print(f"newton_iter diffs:     {len(it_diff)}")
    # iteration delta stats
    if it_diff:
        deltas = [b - a for _, a, b in it_diff]
        print(f"   iter delta: min {min(deltas)}, max {max(deltas)}, "
              f"n_increased {sum(1 for d in deltas if d>0)}, "
              f"n_decreased {sum(1 for d in deltas if d<0)}")
    print(f"converged cases at/over 60-iter cap (new): {len(over60)}")
    for k, it in over60:
        print(f"   {k}  it={it}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", default=None)
    ap.add_argument("--worker-out", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--compare", nargs=2, default=None)
    args = ap.parse_args()
    if args.worker:
        run_worker(args.worker, Path(args.worker_out))
    elif args.compare:
        compare(args.compare[0], args.compare[1])
    elif args.out:
        run_sweep(args.out)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
