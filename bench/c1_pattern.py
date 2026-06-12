#!/usr/bin/env python3
"""C.1 pattern-churn count + extended C.2 bit-identity on the golden case and the
slow-converging sweep cases. Each case runs in a fresh subprocess (so the static
pattern cache starts cold and per-solve counts are clean), with GFOIL_CPAT=1 and
GFOIL_CVERIFY=1. Reports pattern-changes / total-iterations and any memcmp_ok=0.

Run: python3 bench/c1_pattern.py
"""
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FOIL_DIR = REPO / "Smoothed_TEfixed_linear"

# (label, foil_or_None, alpha, ncrit). foil None => golden tests/input.json.
CASES = [
    ("golden NACA0012 a2 nc5", None, 2.0, 5.0),
    ("S4083 (8%) +6 nc9", "S4083 (8%).dat", 6.0, 9.0),
    ("GOE 328 +8 nc9", "GOE 328 AIRFOIL.dat", 8.0, 9.0),
    ("GOE 303 +0 nc9", "GOE 303 (FRIEDRICHSHAFEN G03) AIRFOIL.dat", 0.0, 9.0),
    ("WORTMANN FX 63-137 +6 nc9", "WORTMANN FX 63-137 AIRFOIL.dat", 6.0, 9.0),
    ("AH 79-100 C +6 nc9", "AH 79-100 C AIRFOIL.dat", 6.0, 9.0),
]

WORKER = r'''
import json, sys
sys.path.insert(0, %r)
import GFoil.gfoil_cpp as cpp
foil = %r
alpha = %r
ncrit = %r
if foil is None:
    inp = json.loads(open(%r).read())
    inp["alpha_degrees"] = alpha; inp["ncrit"] = ncrit
else:
    import numpy as np
    from GFoil.inputs import Aerofoil, OperatingConds, Acoustics
    from GFoil.gfoil import _build_input_dict
    xs, ys = [], []
    for line in open(%r, errors="replace"):
        p = line.split()
        if len(p) != 2: continue
        try: x, y = float(p[0]), float(p[1])
        except ValueError: continue
        xs.append(x); ys.append(y)
    foilobj = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0,3.0,0.5]]))
    op = OperatingConds(alpha=alpha, Re=2e6, nCrit=ncrit, rtol=1e-6)
    inp = _build_input_dict(foilobj, op, ac)
r = cpp.run_forward(inp)
print("CONV", r.get("conv",0), r.get("newton_iterations",-1))
'''


def main():
    print(f"{'case':<32s} {'conv':>4s} {'solves':>7s} {'changes':>8s} "
          f"{'scatter':>8s} {'dup_iters':>9s} {'memcmp_fail':>11s}")
    for label, foil, alpha, ncrit in CASES:
        foilpath = str(FOIL_DIR / foil) if foil else ""
        code = WORKER % (str(REPO), foil, alpha, ncrit,
                         str(REPO / "tests" / "input.json"), foilpath)
        env = dict(os.environ, GFOIL_CPAT="1", GFOIL_CVERIFY="1")
        p = subprocess.run([sys.executable, "-c", code],
                           capture_output=True, text=True, env=env)
        cpat = [l for l in p.stderr.splitlines() if l.startswith("[CPAT]")]
        cver = [l for l in p.stderr.splitlines() if l.startswith("[CVERIFY]")]
        total = len(cpat)
        changes = sum(1 for l in cpat if "changed=1" in l)
        scatter = sum(1 for l in cpat if "changed=0" in l)
        memcmp_fail = sum(1 for l in cver if "memcmp_ok=0" in l)
        # duplicate iterations: nonZeros < nnz in any CVERIFY line
        dup = 0
        for l in cver:
            m = re.search(r"nnz=(\d+).*nzA=(\d+)", l)
            if m and int(m.group(2)) < int(m.group(1)):
                dup += 1
        conv = "?"
        for l in p.stdout.splitlines():
            if l.startswith("CONV"):
                conv = l.split()[1]
        print(f"{label:<32s} {conv:>4s} {total:>7d} {changes:>8d} "
              f"{scatter:>8d} {dup:>9d} {memcmp_fail:>11d}")


if __name__ == "__main__":
    main()
