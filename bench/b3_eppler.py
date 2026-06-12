#!/usr/bin/env python3
"""B.3(ii) — EPPLER 399 warm-restart table under each Part-B fix mode.

Drives _call_forward directly (NOT standard_run, so the Part-A driver guard is
not involved — this exercises the C++ solve_coupled fix selected by GFOIL_BFIX).

  off : original criterion  -> warm 5.0<-5.5 instant-accepts (it=0, donor CL)
  A   : honest max(BL,ue)   -> warm re-solves to cold-truth CL (~1.27228), it>0
  B   : skip iter-0 on warm -> warm re-solves to cold-truth CL (~1.27228), it>0

Run: GFOIL_BFIX=<off|A|B> python3 bench/b3_eppler.py
"""
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from GFoil.inputs import Aerofoil, OperatingConds, Acoustics
from GFoil.gfoil import _build_input_dict, _call_forward

FOIL = REPO / "Smoothed_TEfixed_linear" / "EPPLER 399 AIRFOIL.dat"


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


def solve(foil, ac, alpha, donor=None):
    op = OperatingConds(alpha=alpha, Re=2e6, nCrit=5.0, rtol=1e-6)
    inp = _build_input_dict(foil, op, ac, alphaDeg=alpha)
    return _call_forward(inp, prev_result=donor)


def main():
    mode = os.environ.get("GFOIL_BFIX", "B")
    xs, ys = load_dat(FOIL)
    foil = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0, 3.0, 0.5]]))

    print(f"=== GFOIL_BFIX={mode} ===")
    cold50 = solve(foil, ac, 5.0)
    cold55 = solve(foil, ac, 5.5)
    warm = solve(foil, ac, 5.0, donor=cold55)
    print(f"cold 5.0:        it={cold50.newton_iterations:<3d} CL={cold50.CL:.5f}")
    print(f"cold 5.5:        it={cold55.newton_iterations:<3d} CL={cold55.CL:.5f}")
    print(f"warm 5.0 <- 5.5: it={warm.newton_iterations:<3d} CL={warm.CL:.5f}", end="")
    err = abs(warm.CL - 1.27228)
    ok = err < 1e-4
    print(f"  |CL-1.27228|={err:.2e} {'OK(truth)' if ok else 'STALE' if abs(warm.CL-1.32060)<1e-3 else '?'}")


if __name__ == "__main__":
    main()
