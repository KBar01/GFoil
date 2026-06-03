#!/usr/bin/env python3
"""
Boundary extension for the ncrithyst study.

The base sweep (ncrithyst<=0.8) put the convergence optimum AT the grid boundary
(0.8), with the trend still rising. To answer "is more always better, or does it
plateau / hurt?", run the SAME 80-foil grid for higher ncrithyst values, append to
the raw CSV, then re-analyse over the full ordered set.
"""
import sys
import multiprocessing as mp
from pathlib import Path
import pandas as pd

import ncrithyst_study as S   # reuse foil selection, run_one_foil, analysis

EXT_HYST = [1.2, 1.6, 2.0, 3.0]
FULL_HYST = [0.0, 0.1, 0.2, 0.4, 0.8] + EXT_HYST


def main():
    raw = S.OUT_DIR / "ncrithyst_raw.csv"
    base = pd.read_csv(raw)

    # Run only the new ncrithyst values on the identical grid/foils.
    S.HYST = EXT_HYST                      # inherited by forked workers
    foils = S.select_foils()
    print(f"Extension: {len(foils)} foils x ncrithyst{EXT_HYST} "
          f"= {len(foils)*len(S.RES)*len(S.ALPHAS)*len(S.NCRITS)*len(EXT_HYST)} runs")
    tasks = list(enumerate(foils))
    rows = []
    nproc = min(mp.cpu_count(), 10)
    with mp.Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(S.run_one_foil, tasks), 1):
            rows.extend(r)
            print(f"  foil {i}/{len(foils)} done", flush=True)
    ext = pd.DataFrame(rows)

    full = pd.concat([base, ext], ignore_index=True)
    full.to_csv(raw, index=False)
    print(f"Combined raw -> {raw} ({len(full)} rows)")

    # Re-analyse / re-plot / re-write findings over the full ordered ncrithyst set.
    S.HYST = FULL_HYST
    stats = S.analyze(full)
    S.make_plots(stats)
    S.write_findings(stats)
    print("Extension analysis done.")


if __name__ == "__main__":
    main()
