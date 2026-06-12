#!/usr/bin/env python3
"""Phase-1/gate probe for the warm-restart stag seed. cold 6.5 (donor) then
warm 6.5<-6.5 and warm 5.0<-5.5, capturing the run_forward stag trace + the
solve_coupled ENTRY line (GFOIL_DEBUG). Honors GFOIL_NOSTAGSEED / GFOIL_BFIX
from the environment.

Run: GFOIL_DEBUG=1 [GFOIL_NOSTAGSEED=1] python3 bench/stagseed_probe.py
"""
import os
import sys
import tempfile
from contextlib import contextmanager
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
        p = line.split()
        if len(p) != 2:
            continue
        try:
            x, y = float(p[0]), float(p[1])
        except ValueError:
            continue
        xs.append(x); ys.append(y)
    return xs, ys


@contextmanager
def cap():
    sys.stderr.flush()
    saved = os.dup(2)
    tmp = tempfile.TemporaryFile(mode="w+")
    os.dup2(tmp.fileno(), 2)
    try:
        yield tmp
    finally:
        sys.stderr.flush(); os.dup2(saved, 2); os.close(saved)


def solve(foil, ac, alpha, donor=None):
    op = OperatingConds(alpha=alpha, Re=2e6, nCrit=5.0, rtol=1e-6)
    inp = _build_input_dict(foil, op, ac, alphaDeg=alpha)
    with cap() as t:
        r = _call_forward(inp, prev_result=donor)
        t.seek(0); dbg = t.read()
    return r, dbg


def show(label, r, dbg):
    print(f"\n##### {label}: it={r.newton_iterations} "
          f"CL={r.CL if r.converged else 'n/a'}")
    for line in dbg.splitlines():
        if "warm:" in line or "ENTRY" in line:
            print("   " + line)


def main():
    mode = "NOSTAGSEED" if os.environ.get("GFOIL_NOSTAGSEED") else "seed(default)"
    bfix = os.environ.get("GFOIL_BFIX", "B")
    print(f"=== stag seed: {mode}, GFOIL_BFIX={bfix} ===")
    xs, ys = load_dat(FOIL)
    foil = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0, 3.0, 0.5]]))

    d65, _ = solve(foil, ac, 6.5)
    w65, dbg = solve(foil, ac, 6.5, donor=d65)
    show("warm 6.5 <- 6.5", w65, dbg)

    d55, _ = solve(foil, ac, 5.5)
    w50, dbg2 = solve(foil, ac, 5.0, donor=d55)
    show("warm 5.0 <- 5.5", w50, dbg2)
    print(f"   warm 5.0<-5.5 CL={w50.CL:.5f}  |CL-1.27228|={abs(w50.CL-1.27228):.2e}")


if __name__ == "__main__":
    main()
