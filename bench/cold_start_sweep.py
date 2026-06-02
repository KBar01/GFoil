#!/usr/bin/env python3
"""
GFoil cold-start convergence sweep — Phase 0 diagnostic harness.

Establishes the baseline scoreboard against which every solver change is
measured.  Calls the in-process pybind11 solver on the *cold-start path only*
(restart=0, no warm-start / no backstepping continuation), i.e. the equivalent
of GFoil.gfoil._call_forward(inp) with prev_result=None.

Records per case: converged, failure_mode, newton_iterations, wall-time,
CL/CD/CM/OASPL.  Results stream to a JSONL file (resumable) and a summary
scoreboard is printed at the end.

Usage:
    python3 bench/cold_start_sweep.py --run         # run sweep + report
    python3 bench/cold_start_sweep.py --report      # re-print scoreboard only
    python3 bench/cold_start_sweep.py --run --fresh # ignore prior results
    python3 bench/cold_start_sweep.py --run --limit 20   # smoke test
"""

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from bench.foil_select import (build_catalog, select_subset, load_coords,  # noqa: E402
                               geometry, classify, CLASSES, FOIL_DIR)
from GFoil.inputs import Aerofoil, Acoustics, OperatingConds        # noqa: E402
from GFoil.gfoil import _build_input_dict, _call_forward            # noqa: E402

RESULTS_DIR = Path(__file__).resolve().parent / "results"
JSONL_PATH = RESULTS_DIR / "baseline.jsonl"
SUBSET_PATH = RESULTS_DIR / "subset.json"

OBSERVER = np.array([[1.5, 0.0, 1.0]])   # single far-field observer

# Operating-point grid
ALPHAS = [-6.0, -3.0, 0.0, 3.0, 6.0, 10.0]
RES = [2.0e5, 1.0e6, 5.0e6]
NCRITS = [5.0, 9.0, 12.0]
# One forced-transition case per foil (x/c = 0.3 both surfaces).
FORCED_TRANS = (0.3, 0.3)
FORCED_OP = dict(alpha=3.0, Re=1.0e6, nCrit=9.0)


# --------------------------------------------------------------------------- #
# Anchor cases — documented Known Limitations (CLAUDE.md / CHANGELOG.md).      #
# These are fixed regression targets for "did robustness improve".            #
# --------------------------------------------------------------------------- #
def _naca0012_coords():
    """Canonical NACA 0012 used by the golden regression case."""
    inp = json.loads((REPO_ROOT / "tests" / "input.json").read_text())
    return np.asarray(inp["xcoords"]), np.asarray(inp["ycoords"])


def anchor_cases():
    anchors = []
    x12, y12 = _naca0012_coords()
    for a in (2.5, -2.5, 4.7):
        anchors.append(dict(
            anchor=True, foil="NACA 0012 (golden)", cls="thick_symmetric",
            xcoords=x12, ycoords=y12,
            alpha=a, Re=2.0e6, nCrit=5.0, trans=(1.0, 1.0),
            note="nan_lock / ctau-saturation (CLAUDE.md)"))

    x08, y08 = load_coords(REPO_ROOT / "NACA 0008-34.dat")
    anchors.append(dict(
        anchor=True, foil="NACA 0008-34", cls="thick_symmetric",
        xcoords=x08, ycoords=y08,
        alpha=-2.6, Re=2.0e6, nCrit=9.0, trans=(1.0, 1.0),
        note="multi-node BL attractor, period-14 (permanent)"))

    x737, y737 = load_coords(FOIL_DIR / "BOEING 737 MIDSPAN AIRFOIL.dat")
    for a in (-3.1, -3.2):
        anchors.append(dict(
            anchor=True, foil="BOEING 737 MIDSPAN", cls="cambered_gp",
            xcoords=x737, ycoords=y737,
            alpha=a, Re=2.0e6, nCrit=9.0, trans=(1.0, 1.0),
            note="cold-start period-2 oscillation"))
    return anchors


# --------------------------------------------------------------------------- #
# Case construction                                                            #
# --------------------------------------------------------------------------- #
def grid_cases(subset):
    """Full op-grid + one forced-transition case per selected foil."""
    cases = []
    for e in subset:
        x, z = load_coords(Path(e["path"]))
        for a in ALPHAS:
            for re in RES:
                for nc in NCRITS:
                    cases.append(dict(
                        anchor=False, foil=e["name"], cls=e["class"],
                        xcoords=x, ycoords=z,
                        alpha=a, Re=re, nCrit=nc, trans=(1.0, 1.0)))
        cases.append(dict(
            anchor=False, foil=e["name"], cls=e["class"],
            xcoords=x, ycoords=z,
            alpha=FORCED_OP["alpha"], Re=FORCED_OP["Re"], nCrit=FORCED_OP["nCrit"],
            trans=FORCED_TRANS))
    return cases


def case_id(c):
    forced = (c["trans"][0] != 1.0 or c["trans"][1] != 1.0)
    tag = "ANCHOR" if c["anchor"] else ("FT" if forced else "G")
    return (f"{tag}|{c['foil']}|a={c['alpha']:+.2f}|Re={c['Re']:.0e}|"
            f"nc={c['nCrit']:.0f}|tr={c['trans'][0]:.2f}")


# --------------------------------------------------------------------------- #
# Execution                                                                    #
# --------------------------------------------------------------------------- #
def run_case(c):
    """Cold-start solve (restart=0, no warm start).  Returns a result dict."""
    aero = Aerofoil(xcoords=np.asarray(c["xcoords"]).copy(),
                    ycoords=np.asarray(c["ycoords"]).copy())
    op = OperatingConds(alpha=c["alpha"], Re=c["Re"], nCrit=c["nCrit"],
                        ncrithyst=0.0,
                        transition=np.array(c["trans"], dtype=float))
    ac = Acoustics(observerXYZ=OBSERVER)
    inp = _build_input_dict(aero, op, ac, fromRestart=0)

    t0 = time.perf_counter()
    err = ""
    try:
        r = _call_forward(inp)          # prev_result=None -> cold start
        conv = bool(r.converged)
        fm = r.failure_mode
        iters = int(r.newton_iterations)
        cl, cd, cm, oaspl = r.CL, r.CD, r.CM, r.OASPL
    except Exception as e:               # pragma: no cover - defensive
        conv, fm, iters = False, "exception", 0
        cl = cd = cm = oaspl = float("nan")
        err = f"{type(e).__name__}: {e}"
    dt = time.perf_counter() - t0

    return dict(
        case_id=case_id(c), foil=c["foil"], cls=c["cls"], anchor=c["anchor"],
        alpha=c["alpha"], Re=c["Re"], nCrit=c["nCrit"],
        forced=bool(c["trans"][0] != 1.0 or c["trans"][1] != 1.0),
        converged=conv, failure_mode=fm if not conv else "",
        newton_iterations=iters, time_s=dt,
        CL=cl, CD=cd, CM=cm, OASPL=oaspl, error=err,
        note=c.get("note", ""))


def run_sweep(limit=None, fresh=False):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    catalog = build_catalog()
    subset = select_subset(catalog)
    SUBSET_PATH.write_text(json.dumps(
        [{k: e[k] for k in ("name", "class", "t_max", "camber_max", "reflex")}
         for e in subset], indent=2))

    cases = anchor_cases() + grid_cases(subset)
    if limit:
        cases = cases[:limit]

    done = set()
    if JSONL_PATH.exists() and not fresh:
        for line in JSONL_PATH.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["case_id"])
        print(f"Resuming: {len(done)} cases already recorded.")
    elif fresh and JSONL_PATH.exists():
        JSONL_PATH.unlink()

    todo = [c for c in cases if case_id(c) not in done]
    print(f"Total cases: {len(cases)}  |  to run: {len(todo)}\n")

    t_start = time.perf_counter()
    with open(JSONL_PATH, "a") as fh:
        for i, c in enumerate(todo, 1):
            res = run_case(c)
            fh.write(json.dumps(res) + "\n")
            fh.flush()
            flag = "OK " if res["converged"] else "XX "
            if i % 25 == 0 or not res["converged"] or res["anchor"]:
                print(f"  [{i:4d}/{len(todo)}] {flag}{res['case_id']:<60s} "
                      f"it={res['newton_iterations']:3d} "
                      f"{res['time_s']*1e3:6.0f}ms {res['failure_mode']}")
    print(f"\nSweep wall-time: {time.perf_counter() - t_start:.1f}s")


# --------------------------------------------------------------------------- #
# Reporting                                                                    #
# --------------------------------------------------------------------------- #
def load_results():
    if not JSONL_PATH.exists():
        print(f"No results at {JSONL_PATH}. Run with --run first.")
        sys.exit(1)
    return [json.loads(l) for l in JSONL_PATH.read_text().splitlines() if l.strip()]


def _pct(n, d):
    return f"{100.0 * n / d:5.1f}%" if d else "  n/a"


def report():
    rows = load_results()
    grid = [r for r in rows if not r["anchor"]]
    anchors = [r for r in rows if r["anchor"]]

    print("=" * 78)
    print("GFoil COLD-START CONVERGENCE SCOREBOARD")
    print("=" * 78)

    nconv = sum(r["converged"] for r in grid)
    print(f"\nGrid cases: {len(grid)}   cold-converged: {nconv}  "
          f"({_pct(nconv, len(grid))})")

    # Failure-mode breakdown
    print("\nFailure-mode breakdown (non-converged grid cases):")
    modes = {}
    for r in grid:
        if not r["converged"]:
            modes[r["failure_mode"] or "(blank)"] = modes.get(r["failure_mode"] or "(blank)", 0) + 1
    if not modes:
        print("  (none)")
    for m, n in sorted(modes.items(), key=lambda kv: -kv[1]):
        print(f"  {m:<34s} {n:4d}  ({_pct(n, len(grid))})")

    # By geometry class
    print("\nBy geometry class:")
    for c in CLASSES:
        sub = [r for r in grid if r["cls"] == c]
        nc = sum(r["converged"] for r in sub)
        print(f"  {c:<18s} {nc:4d}/{len(sub):<4d} {_pct(nc, len(sub))}")

    # By operating-grid cell
    print("\nBy alpha:")
    for a in ALPHAS:
        sub = [r for r in grid if abs(r["alpha"] - a) < 1e-6 and not r["forced"]]
        nc = sum(r["converged"] for r in sub)
        print(f"  alpha={a:+5.1f}  {nc:4d}/{len(sub):<4d} {_pct(nc, len(sub))}")
    print("\nBy Re:")
    for re in RES:
        sub = [r for r in grid if abs(r["Re"] - re) / re < 1e-6 and not r["forced"]]
        nc = sum(r["converged"] for r in sub)
        print(f"  Re={re:8.0e}  {nc:4d}/{len(sub):<4d} {_pct(nc, len(sub))}")
    print("\nBy nCrit:")
    for ncr in NCRITS:
        sub = [r for r in grid if abs(r["nCrit"] - ncr) < 1e-6 and not r["forced"]]
        nc = sum(r["converged"] for r in sub)
        print(f"  nCrit={ncr:4.0f}  {nc:4d}/{len(sub):<4d} {_pct(nc, len(sub))}")

    forced = [r for r in grid if r["forced"]]
    nf = sum(r["converged"] for r in forced)
    print(f"\nForced-transition cases: {nf}/{len(forced)} {_pct(nf, len(forced))}")

    # Iteration / time stats (converged only)
    conv_iters = [r["newton_iterations"] for r in grid if r["converged"]]
    conv_time = [r["time_s"] for r in grid if r["converged"]]
    if conv_iters:
        print(f"\nConverged Newton iters: mean={statistics.mean(conv_iters):.1f} "
              f"median={statistics.median(conv_iters):.0f} "
              f"max={max(conv_iters)}")
        print(f"Converged wall-time:    mean={statistics.mean(conv_time)*1e3:.0f}ms "
              f"median={statistics.median(conv_time)*1e3:.0f}ms "
              f"total={sum(r['time_s'] for r in grid):.1f}s")

    # Anchor table
    print("\n" + "-" * 78)
    print("DOCUMENTED FAILURE ANCHORS (cold-start):")
    print("-" * 78)
    for r in anchors:
        status = "CONVERGED" if r["converged"] else f"FAILED [{r['failure_mode']}]"
        print(f"  {r['foil']:<22s} a={r['alpha']:+5.2f} nc={r['nCrit']:.0f}  "
              f"{status:<32s} it={r['newton_iterations']:3d}  ({r['note']})")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="store_true", help="run the sweep")
    ap.add_argument("--report", action="store_true", help="print scoreboard")
    ap.add_argument("--fresh", action="store_true", help="discard prior results")
    ap.add_argument("--limit", type=int, default=None, help="cap number of cases")
    args = ap.parse_args()

    if not (args.run or args.report):
        ap.print_help()
        sys.exit(1)
    if args.run:
        run_sweep(limit=args.limit, fresh=args.fresh)
    report()


if __name__ == "__main__":
    main()
