#!/usr/bin/env python3
"""A.3 — re-measure the fwd_run rescue outcome for the 133 cold-failure cases
with the Part-A warm-restart stale-accept guard active.

Only the rescue path is re-run (the 133 cold failures, conv==0 in
cold_start_sweep.csv); the 707 cold successes are untouched. For each case the
full fwd_run continuation runs in its own subprocess (crash isolation) and we
record:
  - rescue_conv_guarded : fwd_run.converged with the guard active
  - guard_fired         : whether "[gfoil] warm-restart stale accept rejected"
                          appeared during the rescue

Outputs:
  - adds a rescue_conv_guarded column to cold_start_sweep.csv (does NOT overwrite
    rescue_conv)
  - appends an addendum section to REPORT.md

Run: source /home/pa20830/NewGradientVal/venv/bin/activate
     python3 bench/rescue_guard.py
     python3 bench/rescue_guard.py --worker "<foil>" <alpha> <ncrit>   # internal
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FOIL_DIR = REPO / "Smoothed_TEfixed_linear"
CSV_PATH = REPO / "bench" / "results" / "cold_start_sweep.csv"
REPORT_PATH = REPO / "bench" / "results" / "REPORT.md"
RE = 2e6
RTOL = 1e-6


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


def worker(foil_name, alpha, ncrit):
    sys.path.insert(0, str(REPO))
    import numpy as np
    from GFoil.inputs import Aerofoil, OperatingConds, Acoustics
    from GFoil.gfoil import fwd_run
    xs, ys = load_dat(FOIL_DIR / foil_name)
    foil = Aerofoil(xcoords=np.array(xs), ycoords=np.array(ys))
    ac = Acoustics(observerXYZ=np.array([[0.0, 3.0, 0.5]]))
    op = OperatingConds(alpha=alpha, Re=RE, nCrit=ncrit, rtol=RTOL)
    r = fwd_run(foil, op, ac)
    print(f"RESULT conv={int(r.converged)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", default=None)
    ap.add_argument("rest", nargs="*")
    args = ap.parse_args()
    if args.worker is not None:
        worker(args.worker, float(args.rest[0]), float(args.rest[1]))
        return

    rows = list(csv.DictReader(open(CSV_PATH)))
    fails = [r for r in rows if r["conv"] == "0"]
    print(f"Re-measuring {len(fails)} cold-failure rescues with guard active ...\n")

    guarded = {}     # (foil, alpha, ncrit) -> (conv, guard_fired)
    for i, r in enumerate(fails, 1):
        key = (r["foil"], r["alpha"], r["ncrit"])
        p = subprocess.run(
            [sys.executable, __file__, "--worker", r["foil"], r["alpha"], r["ncrit"]],
            capture_output=True, text=True)
        conv = 0
        for line in p.stdout.splitlines():
            if line.startswith("RESULT conv="):
                conv = int(line.split("=")[1])
        fired = "warm-restart stale accept rejected" in p.stdout
        guarded[key] = (conv, fired)
        old = r["rescue_conv"]
        flag = ""
        if old == "1" and conv == 0:
            flag = "  <- LOST (was stale accept)"
        elif old != "1" and conv == 1:
            flag = "  <- NEW rescue"
        print(f"[{i:3d}/{len(fails)}] {r['foil']:<42s} a={float(r['alpha']):+.0f} "
              f"nc={float(r['ncrit']):g}  old={old or '0'} new={conv} "
              f"fired={int(fired)}{flag}")

    # --- write CSV with new column -------------------------------------------
    fieldnames = list(rows[0].keys())
    if "rescue_conv_guarded" not in fieldnames:
        fieldnames.append("rescue_conv_guarded")
    for r in rows:
        key = (r["foil"], r["alpha"], r["ncrit"])
        if key in guarded:
            r["rescue_conv_guarded"] = str(guarded[key][0])
        else:
            r["rescue_conv_guarded"] = ""   # cold successes: rescue path not run
    with open(CSV_PATH, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    # --- summary --------------------------------------------------------------
    old_resc = sum(1 for r in fails if r["rescue_conv"] == "1")
    new_resc = sum(1 for k, v in guarded.items() if v[0] == 1)
    lost = [r for r in fails
            if r["rescue_conv"] == "1"
            and guarded[(r["foil"], r["alpha"], r["ncrit"])][0] == 0]
    gained = [r for r in fails
              if r["rescue_conv"] != "1"
              and guarded[(r["foil"], r["alpha"], r["ncrit"])][0] == 1]
    fired_n = sum(1 for v in guarded.values() if v[1])

    print(f"\nprior rescues:   {old_resc}/{len(fails)}")
    print(f"guarded rescues: {new_resc}/{len(fails)}")
    print(f"lost (stale):    {len(lost)}")
    print(f"gained:          {len(gained)}")
    print(f"guard fired in:  {fired_n} rescue runs")

    # --- append addendum to REPORT.md ----------------------------------------
    L = []
    L.append("\n---\n")
    L.append("# Addendum (Part A): rescue re-measurement with warm-restart guard\n")
    L.append("`GFoil/gfoil.py` now rejects a warm-restart iteration-0 accept at a "
             "changed alpha (`_is_stale_warm_accept`) and cold-retries the same "
             "alpha. Re-running the full fwd_run continuation on the 133 cold "
             "failures (`bench/rescue_guard.py`, guard active):\n")
    L.append(f"- prior rescues (unguarded `rescue_conv`): **{old_resc}/{len(fails)}**")
    L.append(f"- corrected rescues (`rescue_conv_guarded`): **{new_resc}/{len(fails)}**")
    L.append(f"- prior \"rescues\" that were stale accepts (lost): **{len(lost)}**")
    L.append(f"- newly-rescued under the guard: **{len(gained)}**")
    L.append(f"- guard fired during **{fired_n}** of the {len(fails)} rescue runs\n")
    L.append("The corrected count replaces the suspect 97. `rescue_conv` is kept "
             "as-is in the CSV; `rescue_conv_guarded` is the new column.\n")

    if lost:
        L.append("## Cases that flipped 1 -> 0 (prior rescue was a stale accept)\n")
        for r in sorted(lost, key=lambda r: (r["foil"], float(r["alpha"]))):
            L.append(f"- `{r['foil']}` alpha={float(r['alpha']):+.0f} "
                     f"nCrit={float(r['ncrit']):g} — {r['failure_mode']}")
        L.append("")
    if gained:
        L.append("## Cases that flipped 0 -> 1 under the guard\n")
        for r in sorted(gained, key=lambda r: (r["foil"], float(r["alpha"]))):
            L.append(f"- `{r['foil']}` alpha={float(r['alpha']):+.0f} "
                     f"nCrit={float(r['ncrit']):g} — {r['failure_mode']}")
        L.append("")

    L.append("## Full per-case rescue table (cold failures)\n")
    L.append("| foil | alpha | nCrit | mode | rescue_conv | rescue_conv_guarded | guard_fired |")
    L.append("|---|---|---|---|---|---|---|")
    for r in sorted(fails, key=lambda r: (r["failure_mode"], r["foil"], float(r["alpha"]))):
        key = (r["foil"], r["alpha"], r["ncrit"])
        conv, fired = guarded[key]
        L.append(f"| `{r['foil']}` | {float(r['alpha']):+.0f} | {float(r['ncrit']):g} "
                 f"| {r['failure_mode']} | {r['rescue_conv'] or '0'} | {conv} | {int(fired)} |")
    L.append("")

    with open(REPORT_PATH, "a") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\naddendum appended to {REPORT_PATH}")
    print(f"rescue_conv_guarded column written to {CSV_PATH}")


if __name__ == "__main__":
    main()
