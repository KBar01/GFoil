#!/usr/bin/env python3
"""Gate (3) — re-measure the 133-case rescue with the warm-restart stag seed,
capturing rescue outcome AND per-rescue Newton-iteration cost. Honors
GFOIL_NOSTAGSEED (pass --noseed) so before/after run on the same binary.

Per case (own subprocess, crash isolation): monkeypatches gfoil_cpp.run_forward
to sum the converged-call Newton iterations and count solves across the whole
fwd_run continuation, then records (converged, conv_iters, n_solves).

Usage:
  python3 bench/rescue_cost.py --out <csv> [--noseed]
  python3 bench/rescue_cost.py --compare <before.csv> <after.csv>
  python3 bench/rescue_cost.py --worker "<foil>" <alpha> <ncrit> [--noseed]  # internal
"""
import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FOIL_DIR = REPO / "Smoothed_TEfixed_linear"
CSV_PATH = REPO / "bench" / "results" / "cold_start_sweep.csv"
RE, RTOL = 2e6, 1e-6


def load_dat(path):
    xs, ys = [], []
    for line in path.read_text(errors="replace").splitlines():
        p = line.split()
        if len(p) != 2:
            continue
        try:
            x, y = float(p[0]), float(p[1])
        except ValueError:
            continue
        xs.append(x); ys.append(y)
    return xs, ys


def worker(foil_name, alpha, ncrit):
    sys.path.insert(0, str(REPO))
    import numpy as np
    import GFoil.gfoil_cpp as cpp
    acc = {"it": 0, "n": 0}
    _orig = cpp.run_forward
    def wrapped(inp, *a):
        r = _orig(inp, *a)
        acc["n"] += 1
        if r.get("conv", 0) == 1:
            acc["it"] += int(r.get("newton_iterations", 0))
        return r
    cpp.run_forward = wrapped
    from GFoil.inputs import Aerofoil, OperatingConds, Acoustics
    from GFoil.gfoil import fwd_run
    xs, ys = load_dat(FOIL_DIR / foil_name)
    foil = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0, 3.0, 0.5]]))
    op = OperatingConds(alpha=alpha, Re=RE, nCrit=ncrit, rtol=RTOL)
    r = fwd_run(foil, op, ac)
    print(f"RESULT conv={int(r.converged)} conv_iters={acc['it']} n_solves={acc['n']}")


def run(out_csv, noseed):
    rows = list(csv.DictReader(open(CSV_PATH)))
    fails = [r for r in rows if r["conv"] == "0"]
    print(f"Re-measuring {len(fails)} rescues "
          f"({'NOSTAGSEED' if noseed else 'seed'}) ...")
    env = dict(os.environ)
    if noseed:
        env["GFOIL_NOSTAGSEED"] = "1"
    else:
        env.pop("GFOIL_NOSTAGSEED", None)
    out = []
    for i, r in enumerate(fails, 1):
        cmd = [sys.executable, __file__, "--worker", r["foil"], r["alpha"], r["ncrit"]]
        p = subprocess.run(cmd, capture_output=True, text=True, env=env)
        conv, it, ns = 0, 0, 0
        for line in p.stdout.splitlines():
            if line.startswith("RESULT"):
                parts = dict(kv.split("=") for kv in line.split()[1:])
                conv, it, ns = int(parts["conv"]), int(parts["conv_iters"]), int(parts["n_solves"])
        out.append({"foil": r["foil"], "alpha": r["alpha"], "ncrit": r["ncrit"],
                    "conv": conv, "conv_iters": it, "n_solves": ns})
        if i % 20 == 0:
            print(f"  {i}/{len(fails)}")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["foil", "alpha", "ncrit", "conv", "conv_iters", "n_solves"])
        w.writeheader(); w.writerows(out)
    resc = sum(r["conv"] for r in out)
    tot = sum(r["conv_iters"] for r in out if r["conv"])
    print(f"rescued {resc}/{len(out)}; total converged-iters over rescued cases: {tot} -> {out_csv}")


def compare(before, after):
    b = {(r["foil"], r["alpha"], r["ncrit"]): r for r in csv.DictReader(open(before))}
    a = {(r["foil"], r["alpha"], r["ncrit"]): r for r in csv.DictReader(open(after))}
    keys = sorted(set(b) & set(a))
    rb = sum(int(b[k]["conv"]) for k in keys)
    ra = sum(int(a[k]["conv"]) for k in keys)
    # iteration totals over cases rescued in BOTH (fair comparison)
    both = [k for k in keys if int(b[k]["conv"]) and int(a[k]["conv"])]
    tb = sum(int(b[k]["conv_iters"]) for k in both)
    ta = sum(int(a[k]["conv_iters"]) for k in both)
    print(f"=== rescue: before {rb}/{len(keys)}  after {ra}/{len(keys)} ===")
    flips = [(k, int(b[k]["conv"]), int(a[k]["conv"])) for k in keys
             if int(b[k]["conv"]) != int(a[k]["conv"])]
    for k, bc, ac_ in flips:
        print(f"   flip {k}: {bc} -> {ac_}")
    print(f"=== converged-iter cost over {len(both)} commonly-rescued cases: "
          f"before {tb}  after {ta}  ({100.0*(ta-tb)/tb:+.1f}%) ===")
    improved = [(k, int(b[k]["conv_iters"]), int(a[k]["conv_iters"])) for k in both
                if int(a[k]["conv_iters"]) < int(b[k]["conv_iters"])]
    worse = [(k, int(b[k]["conv_iters"]), int(a[k]["conv_iters"])) for k in both
             if int(a[k]["conv_iters"]) > int(b[k]["conv_iters"])]
    print(f"   cases cheaper after: {len(improved)}, costlier after: {len(worse)}")
    for k, bi, ai in sorted(improved, key=lambda t: t[2]-t[1])[:12]:
        print(f"     {k}: {bi} -> {ai} iters")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", default=None)
    ap.add_argument("rest", nargs="*")
    ap.add_argument("--out", default=None)
    ap.add_argument("--noseed", action="store_true")
    ap.add_argument("--compare", nargs=2, default=None)
    args = ap.parse_args()
    if args.worker is not None:
        worker(args.worker, float(args.rest[0]), float(args.rest[1]))
    elif args.compare:
        compare(args.compare[0], args.compare[1])
    elif args.out:
        run(args.out, args.noseed)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
