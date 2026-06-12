#!/usr/bin/env python3
"""C.4(iv) — wall-time: scatter (default) vs GFOIL_NOSCATTER (pre-C control), on
the same binary. Golden case (30 reps) + the 10 slowest sweep cases (7 reps).
Each (case, mode) runs in a fresh subprocess so the env takes effect and the
pattern cache warms identically; the steady-state (post-warmup) median is what
the scatter path improves.

Run: python3 bench/c4_timing.py
"""
import os
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FOIL_DIR = REPO / "Smoothed_TEfixed_linear"

# (label, foil_or_None, alpha, ncrit, reps)
CASES = [
    ("golden NACA0012 a2 nc5", None, 2.0, 5.0, 30),
    ("S4083 (8%) +6 nc9", "S4083 (8%).dat", 6.0, 9.0, 7),
    ("BELL-WORTMANN FX69 +4 nc5", "BELL-WORTMANN FX 69-H-083 AIRFOIL.dat", 4.0, 5.0, 7),
    ("GOE 303 +0 nc9", "GOE 303 (FRIEDRICHSHAFEN G03) AIRFOIL.dat", 0.0, 9.0, 7),
    ("GOE 328 +8 nc9", "GOE 328 AIRFOIL.dat", 8.0, 9.0, 7),
    ("GOE 458 -2 nc9", "GOE 458 AIRFOIL.dat", -2.0, 9.0, 7),
    ("NASA MS(1)-0313 +8 nc9", "NASA-LANGLEY MS(1)-0313 AIRFOIL.dat", 8.0, 9.0, 7),
    ("WORTMANN FX 63-137 +6 nc9", "WORTMANN FX 63-137 AIRFOIL.dat", 6.0, 9.0, 7),
    ("AH 79-100 C +6 nc9", "AH 79-100 C AIRFOIL.dat", 6.0, 9.0, 7),
    ("MH 93 15.98% -4 nc5", "MH 93  15.98%.dat", -4.0, 5.0, 7),
    ("E374 -2 nc9", "E374.dat", -2.0, 9.0, 7),
]

WORKER = r'''
import json, sys, time, statistics
sys.path.insert(0, %r)
import GFoil.gfoil_cpp as cpp
foil, alpha, ncrit, reps = %r, %r, %r, %d
if foil is None:
    inp = json.loads(open(%r).read()); inp["alpha_degrees"]=alpha; inp["ncrit"]=ncrit
else:
    import numpy as np
    from GFoil.inputs import Aerofoil, OperatingConds, Acoustics
    from GFoil.gfoil import _build_input_dict
    xs, ys = [], []
    for line in open(%r, errors="replace"):
        p = line.split()
        if len(p)!=2: continue
        try: x,y=float(p[0]),float(p[1])
        except ValueError: continue
        xs.append(x); ys.append(y)
    fo = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0,3.0,0.5]]))
    op = OperatingConds(alpha=alpha, Re=2e6, nCrit=ncrit, rtol=1e-6)
    inp = _build_input_dict(fo, op, ac)
for _ in range(3): cpp.run_forward(inp)        # warmup
ts=[]
for _ in range(reps):
    t0=time.perf_counter(); cpp.run_forward(inp); ts.append((time.perf_counter()-t0)*1e3)
print("MS", statistics.median(ts))
'''


def time_case(foil, alpha, ncrit, reps, noscatter):
    foilpath = str(FOIL_DIR / foil) if foil else ""
    code = WORKER % (str(REPO), foil, alpha, ncrit, reps,
                     str(REPO / "tests" / "input.json"), foilpath)
    env = dict(os.environ)
    if noscatter:
        env["GFOIL_NOSCATTER"] = "1"
    else:
        env.pop("GFOIL_NOSCATTER", None)
    p = subprocess.run([sys.executable, "-c", code],
                       capture_output=True, text=True, env=env)
    for l in p.stdout.splitlines():
        if l.startswith("MS"):
            return float(l.split()[1])
    return float("nan")


def main():
    print(f"{'case':<30s} {'noscatter(ms)':>14s} {'scatter(ms)':>12s} "
          f"{'delta%':>8s}")
    for label, foil, alpha, ncrit, reps in CASES:
        # interleave to reduce drift: ns, sc, ns, sc -> take min of two each
        ns1 = time_case(foil, alpha, ncrit, reps, True)
        sc1 = time_case(foil, alpha, ncrit, reps, False)
        ns2 = time_case(foil, alpha, ncrit, reps, True)
        sc2 = time_case(foil, alpha, ncrit, reps, False)
        ns = min(ns1, ns2); sc = min(sc1, sc2)
        d = 100.0 * (sc - ns) / ns
        print(f"{label:<30s} {ns:>14.1f} {sc:>12.1f} {d:>+7.1f}%")


if __name__ == "__main__":
    main()
