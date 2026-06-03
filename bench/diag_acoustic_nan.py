#!/usr/bin/env python3
"""Phase A0 diagnostic: replay the cold-start acoustic_nan cases and capture the
WPS-intermediate trace (GFOIL_DEBUG) to identify which term first goes
non-finite.  Reads the case list from a prior sweep's JSONL (failure_mode ==
acoustic_nan), or uses a hard-coded fallback list.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from bench.foil_select import load_coords, FOIL_DIR          # noqa: E402
from GFoil.inputs import Aerofoil, OperatingConds, Acoustics  # noqa: E402
from GFoil.gfoil import _build_input_dict                     # noqa: E402
import GFoil.gfoil_cpp as cpp                                  # noqa: E402

OBSERVER = np.array([[1.5, 0.0, 1.0]])


def case_list(tag="final"):
    p = Path(__file__).resolve().parent / "results" / f"sweep_{tag}.jsonl"
    rows = [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
    return [r for r in rows if r.get("failure_mode") == "acoustic_nan"]


def foil_path(name):
    # anchors use root files; otherwise the smoothed folder
    for cand in (REPO / f"{name}.dat", FOIL_DIR / f"{name}.dat",
                 FOIL_DIR / f"{name} AIRFOIL.dat"):
        if cand.exists():
            return cand
    # fuzzy: match stem prefix in the folder
    for f in FOIL_DIR.glob("*.dat"):
        if f.stem == name:
            return f
    raise FileNotFoundError(name)


def run(case, rtol=None, debug=True):
    x, z = load_coords(foil_path(case["foil"]))
    aero = Aerofoil(xcoords=x.copy(), ycoords=z.copy())
    op = OperatingConds(alpha=case["alpha"], Re=case["Re"], nCrit=case["nCrit"])
    ac = Acoustics(observerXYZ=OBSERVER)
    inp = _build_input_dict(aero, op, ac, fromRestart=0)
    if rtol is not None:
        inp["rtol"] = rtol
    r = cpp.run_forward(inp)
    return r


if __name__ == "__main__":
    cases = case_list()
    print(f"{len(cases)} acoustic_nan cases\n", file=sys.stderr)
    os.environ["GFOIL_DEBUG"] = "1"
    summary = []
    for c in cases:
        tag = f"{c['foil']} a={c['alpha']:+.1f} Re={c['Re']:.0e} nc={c['nCrit']:.0f}"
        print(f"\n===== {tag} =====", file=sys.stderr, flush=True)
        r = run(c)
        summary.append((tag, r.get("conv"), r.get("failure_mode", ""),
                        r.get("OASPL", float("nan"))))
    print("\n\n==== SUMMARY ====")
    for tag, conv, fm, oaspl in summary:
        print(f"  conv={conv} fm={fm:<14s} OASPL={oaspl}  {tag}")
