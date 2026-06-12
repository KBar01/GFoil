#!/usr/bin/env python3
"""A.2 — reproduce the EPPLER 399 warm-restart ground-truth table with the
Part-A stale-accept guard active.

Ground truth (REPORT.md "warm-restart stale accept", nCrit=5, Re=2e6, rtol=1e-6):
    cold 5.0          -> it=27, CL=1.27228   (truth)
    cold 5.5          -> it=30, CL=1.32060   (donor)
    warm 5.0 <- 5.5   -> it=0,  CL=1.32056   (the defect: donor echoed)
    warm 6.5 <- 6.5   -> it=8,  CL=1.40848   (same-alpha; must NOT be rejected)

The guard (_is_stale_warm_accept) is the exact predicate standard_run uses in
its forward-stepping loop. This script drives _call_forward directly to show the
raw run_forward behaviour, then applies the guard and verifies the cold-truth CL
is recovered, and that a same-alpha restart is left untouched.

Run: source /home/pa20830/NewGradientVal/venv/bin/activate
     python3 bench/eppler399_guard.py
"""
import os
import sys
from pathlib import Path

import numpy as np

# This script demonstrates the Part-A Python guard, which catches a C++ stale
# accept. With Part-B Fix B active (the default build) the C++ no longer PRODUCES
# a stale accept, so to exercise the guard we reproduce the defect with
# GFOIL_BFIX=off (read per solve_coupled call). Set before importing the module.
os.environ["GFOIL_BFIX"] = "off"

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from GFoil.inputs import Aerofoil, OperatingConds, Acoustics
from GFoil.gfoil import _build_input_dict, _call_forward, _is_stale_warm_accept

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
    xs, ys = load_dat(FOIL)
    foil = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0, 3.0, 0.5]]))

    print("=== EPPLER 399, nCrit=5, Re=2e6, rtol=1e-6 ===\n")

    cold50 = solve(foil, ac, 5.0)
    print(f"cold 5.0:          it={cold50.newton_iterations:<3d} CL={cold50.CL:.5f}  <- truth")
    cold55 = solve(foil, ac, 5.5)
    print(f"cold 5.5:          it={cold55.newton_iterations:<3d} CL={cold55.CL:.5f}  <- donor")

    warm = solve(foil, ac, 5.0, donor=cold55)
    print(f"warm 5.0 <- 5.5:   it={warm.newton_iterations:<3d} CL={warm.CL:.5f}  "
          f"<- raw run_forward (stale accept)")

    # --- guard active ---------------------------------------------------------
    print("\n--- guard active ---")
    stale = _is_stale_warm_accept(warm, 5.0, float(cold55.alpha))
    print(f"_is_stale_warm_accept(warm 5.0 <- 5.5) = {stale}")
    assert stale, "guard FAILED to flag the documented stale accept"
    # Mirror the standard_run print + cold retry.
    print(f"[gfoil] warm-restart stale accept rejected at alpha={5.0:.4f} "
          f"(donor {float(cold55.alpha):.4f}); cold retry")
    guarded = solve(foil, ac, 5.0)  # cold retry
    print(f"warm 5.0 <- 5.5 (guarded): it={guarded.newton_iterations:<3d} "
          f"CL={guarded.CL:.5f}")
    err = abs(guarded.CL - 1.27228)
    print(f"  |CL - 1.27228| = {err:.2e}  ({'PASS' if err < 1e-4 else 'FAIL'} < 1e-4)")
    assert err < 1e-4, "guarded result did not recover cold-truth CL"

    # --- same-alpha must NOT be rejected -------------------------------------
    print("\n--- same-alpha restart (must NOT be rejected) ---")
    cold65 = solve(foil, ac, 6.5)
    warm_same = solve(foil, ac, 6.5, donor=cold65)
    print(f"cold 6.5:          it={cold65.newton_iterations:<3d} CL={cold65.CL:.5f}")
    print(f"warm 6.5 <- 6.5:   it={warm_same.newton_iterations:<3d} CL={warm_same.CL:.5f}")
    same_stale = _is_stale_warm_accept(warm_same, 6.5, float(cold65.alpha))
    print(f"_is_stale_warm_accept(warm 6.5 <- 6.5) = {same_stale}  "
          f"({'PASS' if not same_stale else 'FAIL'}: same-alpha not rejected)")
    assert not same_stale, "guard wrongly rejected a same-alpha restart"

    print("\nA.2 OK: stale accept caught + cold-truth recovered; same-alpha left allowed.")


if __name__ == "__main__":
    main()
