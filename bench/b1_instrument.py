#!/usr/bin/env python3
"""B.1 Phase-1 instrumentation driver (run with GFOIL_DEBUG=1).

Three EPPLER 399 cases (nCrit=5, Re=2e6, rtol=1e-6):
  a. warm 5.0 <- 5.5  : expect BL_rms < rtol, ue_rms >> rtol  (the smoking gun)
  b. warm 6.5 <- 6.5  : same-alpha; explain the it=8 anomaly via the ENTRY (warm)
                        vs CONVERGED (donor) config diff
  c. cold 5.0         : control; both RMS large at entry

C++ writes [GFOIL_DEBUG] lines to fd 2; we capture fd 2 around each call so the
ENTRY/CONVERGED lines are attributed to the right solve.

Run: GFOIL_DEBUG=1 python3 bench/b1_instrument.py
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
        parts = line.split()
        if len(parts) != 2:
            continue
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        xs.append(x); ys.append(y)
    return xs, ys


@contextmanager
def capture_fd2():
    """Capture C-level fd 2 (std::cerr) into a string."""
    sys.stderr.flush()
    saved = os.dup(2)
    tmp = tempfile.TemporaryFile(mode="w+")
    os.dup2(tmp.fileno(), 2)
    try:
        yield tmp
    finally:
        sys.stderr.flush()
        os.dup2(saved, 2)
        os.close(saved)


def solve(foil, ac, alpha, donor=None):
    op = OperatingConds(alpha=alpha, Re=2e6, nCrit=5.0, rtol=1e-6)
    inp = _build_input_dict(foil, op, ac, alphaDeg=alpha)
    with capture_fd2() as tmp:
        r = _call_forward(inp, prev_result=donor)
        tmp.seek(0)
        dbg = tmp.read()
    return r, dbg


def show(label, r, dbg):
    print(f"\n##### {label} : it={r.newton_iterations} "
          f"CL={r.CL if r.converged else 'n/a'}")
    for line in dbg.splitlines():
        if "GFOIL_DEBUG" in line:
            print("   " + line)


def main():
    xs, ys = load_dat(FOIL)
    foil = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0, 3.0, 0.5]]))

    print("=== B.1 case (a): warm 5.0 <- 5.5 (smoking gun) ===")
    donor55, _ = solve(foil, ac, 5.5)
    warm50, dbg = solve(foil, ac, 5.0, donor=donor55)
    show("warm 5.0 <- 5.5", warm50, dbg)

    print("\n\n=== B.1 case (b): warm 6.5 <- 6.5 (same-alpha it=8 anomaly) ===")
    donor65, dbg_d = solve(foil, ac, 6.5)
    show("cold 6.5 (donor)", donor65, dbg_d)
    warm65, dbg_w = solve(foil, ac, 6.5, donor=donor65)
    show("warm 6.5 <- 6.5", warm65, dbg_w)

    print("\n\n=== B.1 case (c): cold 5.0 (control) ===")
    cold50, dbg_c = solve(foil, ac, 5.0)
    show("cold 5.0", cold50, dbg_c)


if __name__ == "__main__":
    main()
