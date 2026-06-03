#!/usr/bin/env python3
"""
ncrithyst convergence study
===========================

Question: does the free-transition hysteresis margin `ncrithyst` (used in exactly
one place, update_transition.cpp ~line 168) measurably help the bare Newton solve
converge, and over what range of values?  The hypothesis is that a non-zero
`ncrithyst` damps single-node transition-front retreat and therefore reduces
"transition_front_oscillation" / "no_convergence" failures.

Controlled comparison: for each (aerofoil, Re, alpha, nCrit) CASE we run the SAME
case once per `ncrithyst` value and compare outcomes.

-----------------------------------------------------------------------------
IMPORTANT — build provenance
-----------------------------------------------------------------------------
The primary working tree (/home/pa20830/GFoil) had `ncrithyst` REMOVED from the
user-facing API earlier (it is now a hard-coded constant and the input-dict key
is gone), so its built `gfoil_cpp` cannot sweep the parameter.  This study
therefore imports the pybind module from a git worktree checked out at HEAD
(commit good_state), where `ncrithyst` is still threaded through
OperatingConds.ncrithyst -> input dict -> gfoil_fwd_bindings.cpp.  NO library
source is modified; update_transition.cpp etc. are exactly the committed code.

-----------------------------------------------------------------------------
Design
-----------------------------------------------------------------------------
* Geometries: a SEEDED (seed=0) random subset of N_FOILS aerofoils from the main
  repo's Smoothed_TEfixed_linear/ (Tecplot .dat), loaded with the repo's
  header-tolerant foil_select.load_coords convention.
* Grid (full factorial):
    ncrithyst in {0.0, 0.1, 0.2, 0.4, 0.8}   (0.0 = feature effectively OFF)
    Re        in {2e5, 5e5, 1e6, 3e6, 6e6, 9e6}
    alpha     in {0, 2, 4, 6, 8, 10, 12} deg
    nCrit     in {6, 9, 12}
    FREE transition only (transition=[1,1]); ncrithyst is a no-op for tripped BLs.
* Fixed: Ma=0, rho=1.225, nu=1.5e-5, chord=1, default panels, single far-field
  observer (convergence only — acoustic values are not the subject).
* Per run we call the BARE single Newton solve (_build_input_dict + cpp.run_forward),
  NOT fwd_run/standard_run, because standard_run performs alpha-backstepping
  continuation that would rescue bare failures and mask ncrithyst's effect.
* Every solve is wrapped in try/except; crashes are logged as failure_mode='exception'.

Outputs (written next to this script):
  ncrithyst_raw.csv          one row per run
  summary_convergence.csv    convergence % by ncrithyst, overall + by Re/alpha band
  failure_tally.csv          failure_mode counts per ncrithyst
  paired_iterations.csv      matched-converged newton_iterations per ncrithyst
  plot_conv_vs_hyst.png      (a) convergence rate vs ncrithyst
  plot_conv_by_Re.png        (b) convergence rate vs ncrithyst, faceted by Re
  plot_failmodes.png         (c) failure_mode stacked bar per ncrithyst
  findings.md                ~1-page written summary
"""

import os
import sys
import glob
import json
import random
import argparse
import traceback
from pathlib import Path

import numpy as np

# --- build provenance: import the HEAD worktree that still exposes ncrithyst ---
WORKTREE = "/home/pa20830/gfoil_ncrithyst_wt"
MAIN_REPO = "/home/pa20830/GFoil"
sys.path.insert(0, WORKTREE)                 # GFoil package WITH ncrithyst
sys.path.insert(0, os.path.join(MAIN_REPO, "bench"))  # foil_select.load_coords

import GFoil.gfoil_cpp as cpp                      # noqa: E402
from GFoil.inputs import Aerofoil, Acoustics, OperatingConds  # noqa: E402
from GFoil.gfoil import _build_input_dict           # noqa: E402
from foil_select import load_coords                 # noqa: E402

OUT_DIR   = Path(__file__).resolve().parent
FOIL_DIR  = Path(MAIN_REPO) / "Smoothed_TEfixed_linear"

# ----------------------------- study grid ---------------------------------- #
SEED       = 0
N_FOILS    = 80
HYST       = [0.0, 0.1, 0.2, 0.4, 0.8]
RES        = [2e5, 5e5, 1e6, 3e6, 6e6, 9e6]
ALPHAS     = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
NCRITS     = [6.0, 9.0, 12.0]
OBSERVER   = np.array([[1.5, 0.0, 1.0]])
MA, RHO, NU, CHORD = 0.0, 1.225, 1.5e-5, 1.0


# --------------------------- foil selection -------------------------------- #
def select_foils():
    """Seeded random subset of valid foils (pass Aerofoil() validation)."""
    files = sorted(glob.glob(str(FOIL_DIR / "*.dat")))
    rng = random.Random(SEED)
    rng.shuffle(files)
    chosen = []
    for f in files:
        try:
            x, z = load_coords(f)
            Aerofoil(xcoords=x.copy(), ycoords=z.copy())  # validate only
            chosen.append(f)
        except Exception:
            continue
        if len(chosen) >= N_FOILS:
            break
    return chosen


# --------------------------- single foil sweep ----------------------------- #
def run_one_foil(args):
    """Run the full (Re, alpha, nCrit, ncrithyst) sub-grid for one foil.

    Returns a list of row dicts. Loads its own coords so workers stay lean.
    """
    foil_idx, path = args
    foil_name = Path(path).stem
    rows = []
    try:
        x, z = load_coords(path)
        foil = Aerofoil(xcoords=x.copy(), ycoords=z.copy(),
                        chord=CHORD, span=2.0)
    except Exception as e:
        return [dict(foil_idx=foil_idx, foil=foil_name, Re=np.nan, alpha=np.nan,
                     nCrit=np.nan, ncrithyst=np.nan, converged=0,
                     failure_mode="foil_load_error", newton_iterations=0,
                     CL=np.nan, CD=np.nan, error=str(e))]

    ac = Acoustics(observerXYZ=OBSERVER)
    for Re in RES:
        for al in ALPHAS:
            for nc in NCRITS:
                for h in HYST:
                    row = dict(foil_idx=foil_idx, foil=foil_name, Re=Re, alpha=al,
                               nCrit=nc, ncrithyst=h, converged=0,
                               failure_mode="", newton_iterations=0,
                               CL=np.nan, CD=np.nan, error="")
                    try:
                        op = OperatingConds(alpha=al, Re=Re, Ma=MA, rho=RHO, nu=NU,
                                            nCrit=nc, ncrithyst=h,
                                            transition=np.array([1.0, 1.0]))
                        inp = _build_input_dict(foil, op, ac, fromRestart=0,
                                                verbose=False)
                        r = cpp.run_forward(inp)   # BARE single Newton solve
                        row["converged"]         = int(r.get("conv", 0))
                        row["failure_mode"]      = r.get("failure_mode", "")
                        row["newton_iterations"] = int(r.get("newton_iterations", 0))
                        if row["converged"]:
                            row["CL"] = float(r.get("CL", np.nan))
                            row["CD"] = float(r.get("CD", np.nan))
                    except Exception as e:
                        row["failure_mode"] = "exception"
                        row["error"] = f"{type(e).__name__}: {e}"
                    rows.append(row)
    return rows


def run_sweep():
    import multiprocessing as mp
    import pandas as pd

    foils = select_foils()
    print(f"Selected {len(foils)} valid foils (seed={SEED}).")
    n_cases = len(foils) * len(RES) * len(ALPHAS) * len(NCRITS)
    print(f"Cases (foil x Re x alpha x nCrit): {n_cases};  "
          f"runs (x{len(HYST)} ncrithyst): {n_cases * len(HYST)}")

    tasks = list(enumerate(foils))
    all_rows = []
    nproc = min(mp.cpu_count(), 10)
    print(f"Running on {nproc} processes ...")
    with mp.Pool(nproc) as pool:
        for i, rows in enumerate(pool.imap_unordered(run_one_foil, tasks), 1):
            all_rows.extend(rows)
            print(f"  foil {i}/{len(foils)} done ({len(rows)} runs)", flush=True)

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_DIR / "ncrithyst_raw.csv", index=False)
    print(f"Wrote {OUT_DIR / 'ncrithyst_raw.csv'}  ({len(df)} rows)")
    return df


# ------------------------------- analysis ---------------------------------- #
def re_band(re):
    if re <= 5e5:
        return "low (<=5e5)"
    if re <= 1e6:
        return "mid (1e6)"
    return "high (>=3e6)"


def alpha_band(a):
    if a <= 2:
        return "low (0-2)"
    if a <= 6:
        return "mid (4-6)"
    return "high (8-12)"


def analyze(df):
    import pandas as pd

    # Derived: aero convergence. failure_mode=='acoustic_nan' means the AERO
    # Newton solve converged but the downstream noise model returned non-finite
    # OASPL (ncrithyst is an aero-only knob), so it counts as an aero success.
    df["aero_conv"] = ((df["converged"] == 1) |
                       (df["failure_mode"] == "acoustic_nan")).astype(int)
    df["case_id"] = (df["foil_idx"].astype(str) + "|" + df["Re"].astype(str) +
                     "|" + df["alpha"].astype(str) + "|" + df["nCrit"].astype(str))
    df["re_band"] = df["Re"].map(re_band)
    df["alpha_band"] = df["alpha"].map(alpha_band)

    # ---- (1) convergence % by ncrithyst, overall + by band -----------------
    rows = []
    for h in HYST:
        sub = df[df["ncrithyst"] == h]
        rows.append(dict(group="overall", value="all", ncrithyst=h,
                         n=len(sub),
                         full_conv_pct=100 * sub["converged"].mean(),
                         aero_conv_pct=100 * sub["aero_conv"].mean()))
    for band in sorted(df["re_band"].unique()):
        for h in HYST:
            sub = df[(df["ncrithyst"] == h) & (df["re_band"] == band)]
            rows.append(dict(group="Re_band", value=band, ncrithyst=h, n=len(sub),
                             full_conv_pct=100 * sub["converged"].mean(),
                             aero_conv_pct=100 * sub["aero_conv"].mean()))
    for band in sorted(df["alpha_band"].unique()):
        for h in HYST:
            sub = df[(df["ncrithyst"] == h) & (df["alpha_band"] == band)]
            rows.append(dict(group="alpha_band", value=band, ncrithyst=h, n=len(sub),
                             full_conv_pct=100 * sub["converged"].mean(),
                             aero_conv_pct=100 * sub["aero_conv"].mean()))
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "summary_convergence.csv", index=False)

    # ---- failure_mode tally per ncrithyst ----------------------------------
    df["fmode"] = df["failure_mode"].replace("", "converged")
    tally = (df.groupby(["ncrithyst", "fmode"]).size()
             .unstack(fill_value=0).reset_index())
    tally.to_csv(OUT_DIR / "failure_tally.csv", index=False)

    # ---- paired analysis: matched-converged subset -------------------------
    # Cases where ALL ncrithyst values aero-converged -> compare iterations.
    piv = df.pivot_table(index="case_id", columns="ncrithyst",
                         values="aero_conv", aggfunc="first")
    matched = piv.index[(piv[HYST] == 1).all(axis=1)]
    it = df[df["case_id"].isin(matched)]
    paired = (it.groupby("ncrithyst")["newton_iterations"]
              .agg(["count", "mean", "median"]).reset_index())
    paired.to_csv(OUT_DIR / "paired_iterations.csv", index=False)

    # ---- rescues / regressions vs ncrithyst=0 ------------------------------
    conv0 = piv[0.0]
    rescue = {}      # h -> # cases conv at h but NOT at 0.0
    regress = {}     # h -> # cases conv at 0.0 but NOT at h
    for h in HYST:
        if h == 0.0:
            continue
        rescue[h] = int(((piv[h] == 1) & (conv0 == 0)).sum())
        regress[h] = int(((piv[h] == 0) & (conv0 == 1)).sum())

    stats = dict(summary=summary, tally=tally, paired=paired,
                 matched_n=len(matched), total_cases=piv.shape[0],
                 rescue=rescue, regress=regress, df=df)
    return stats


# ------------------------------- plots ------------------------------------- #
def make_plots(stats):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = stats["df"]
    summary = stats["summary"]

    # (a) convergence rate vs ncrithyst (overall)
    ov = summary[summary["group"] == "overall"].sort_values("ncrithyst")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ov["ncrithyst"], ov["aero_conv_pct"], "o-", label="aero converged")
    ax.plot(ov["ncrithyst"], ov["full_conv_pct"], "s--", label="full (incl. acoustic)")
    ax.set_xlabel("ncrithyst"); ax.set_ylabel("convergence rate [%]")
    ax.set_title(f"Convergence vs ncrithyst  (n={len(df)//len(HYST)} cases/value)")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig(OUT_DIR / "plot_conv_vs_hyst.png", dpi=120)
    plt.close(fig)

    # (b) faceted by Re
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for re in RES:
        sub = (df[df["Re"] == re].groupby("ncrithyst")["aero_conv"]
               .mean().reset_index())
        ax.plot(sub["ncrithyst"], 100 * sub["aero_conv"], "o-",
                label=f"Re={re:.0e}")
    ax.set_xlabel("ncrithyst"); ax.set_ylabel("aero convergence rate [%]")
    ax.set_title("Convergence vs ncrithyst, by Re")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT_DIR / "plot_conv_by_Re.png", dpi=120)
    plt.close(fig)

    # (c) failure_mode stacked bar per ncrithyst
    tally = stats["tally"].set_index("ncrithyst")
    cols = [c for c in tally.columns if c != "converged"]  # show failures only
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bottom = np.zeros(len(tally))
    for c in cols:
        ax.bar(tally.index.astype(str), tally[c], bottom=bottom, label=c, width=0.6)
        bottom += tally[c].values
    ax.set_xlabel("ncrithyst"); ax.set_ylabel("# failed runs")
    ax.set_title("Failure modes per ncrithyst (failures only)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT_DIR / "plot_failmodes.png", dpi=120)
    plt.close(fig)


# ------------------------------ findings ----------------------------------- #
def write_findings(stats):
    summary = stats["summary"]
    paired = stats["paired"]
    tally = stats["tally"].set_index("ncrithyst")
    ov = summary[summary["group"] == "overall"].sort_values("ncrithyst")

    base = ov[ov["ncrithyst"] == 0.0]["aero_conv_pct"].iloc[0]
    best_row = ov.loc[ov["aero_conv_pct"].idxmax()]
    best_h = best_row["ncrithyst"]
    best_pct = best_row["aero_conv_pct"]
    delta = best_pct - base

    tfo = (tally["transition_front_oscillation"]
           if "transition_front_oscillation" in tally.columns
           else None)

    lines = []
    lines.append("# ncrithyst convergence study — findings\n")
    lines.append(f"*Seed={SEED}, {N_FOILS} foils requested; "
                 f"grid ncrithyst{HYST} x Re{['%.0e'%r for r in RES]} x "
                 f"alpha{ALPHAS} x nCrit{NCRITS}, free transition only. "
                 f"Bare single Newton solve (no backstepping/repanel). "
                 f"Total runs: {len(stats['df'])}.*\n")

    lines.append("## Headline\n")
    verdict = ("**measurably helps**" if delta >= 1.0 else
               "**no meaningful effect (null result)**" if abs(delta) < 1.0 else
               "**hurts**")
    lines.append(f"Across the swept grid, increasing `ncrithyst` from 0.0 changes the "
                 f"aero-convergence rate from **{base:.1f}%** (off) to a maximum of "
                 f"**{best_pct:.1f}%** at ncrithyst={best_h:g} "
                 f"(Δ = {delta:+.1f} percentage points). Conclusion: ncrithyst "
                 f"{verdict} for bare Newton convergence on this sample.\n")

    lines.append("## Convergence rate vs ncrithyst (overall)\n")
    lines.append("| ncrithyst | aero conv % | full conv % | n |")
    lines.append("|---|---|---|---|")
    for _, r in ov.iterrows():
        lines.append(f"| {r['ncrithyst']:g} | {r['aero_conv_pct']:.2f} | "
                     f"{r['full_conv_pct']:.2f} | {int(r['n'])} |")
    lines.append("")

    lines.append("## Rescues and regressions vs ncrithyst=0\n")
    lines.append("Cases (foil,Re,alpha,nCrit) that converge ONLY with ncrithyst>0 "
                 "(rescue) vs ONLY with ncrithyst=0 (regression):\n")
    lines.append("| ncrithyst | rescued (only h>0) | regressed (only h=0) |")
    lines.append("|---|---|---|")
    for h in HYST:
        if h == 0.0:
            continue
        lines.append(f"| {h:g} | {stats['rescue'][h]} | {stats['regress'][h]} |")
    lines.append("")

    lines.append("## transition_front_oscillation — the targeted failure mode\n")
    if tfo is not None:
        lines.append("| ncrithyst | # transition_front_oscillation |")
        lines.append("|---|---|")
        for h in HYST:
            v = int(tfo.get(h, 0)) if h in tfo.index else 0
            lines.append(f"| {h:g} | {v} |")
        t0 = int(tfo.get(0.0, 0))
        lines.append("")
        lines.append(f"At ncrithyst=0 there were **{t0}** transition_front_oscillation "
                     f"failures. This is the exact failure mode the hysteresis margin "
                     f"targets, so its trend is the smoking gun.\n")
    else:
        lines.append("No `transition_front_oscillation` failures occurred anywhere in "
                     "the sweep — the targeted failure mode did not appear on this "
                     "sample, so ncrithyst had nothing to suppress.\n")

    lines.append("## Speed: matched-converged Newton iterations\n")
    lines.append(f"Among the {stats['matched_n']} of {stats['total_cases']} cases that "
                 f"converge at ALL ncrithyst values:\n")
    lines.append("| ncrithyst | mean iters | median iters |")
    lines.append("|---|---|---|")
    for _, r in paired.iterrows():
        lines.append(f"| {r['ncrithyst']:g} | {r['mean']:.2f} | {int(r['median'])} |")
    lines.append("")

    lines.append("## Regime breakdown (aero conv %)\n")
    for grp in ["Re_band", "alpha_band"]:
        g = summary[summary["group"] == grp]
        bands = list(dict.fromkeys(g["value"]))
        lines.append(f"### by {grp}\n")
        lines.append("| " + grp + " | " +
                     " | ".join(f"h={h:g}" for h in HYST) + " |")
        lines.append("|" + "---|" * (len(HYST) + 1))
        for b in bands:
            cells = []
            for h in HYST:
                v = g[(g["value"] == b) & (g["ncrithyst"] == h)]["aero_conv_pct"]
                cells.append(f"{v.iloc[0]:.1f}" if len(v) else "-")
            lines.append(f"| {b} | " + " | ".join(cells) + " |")
        lines.append("")

    lines.append("## Recommendation\n")
    if delta >= 1.0:
        lines.append(f"Recommend defaulting `ncrithyst = {best_h:g}` "
                     f"(best convergence, Δ{delta:+.1f} pts vs off). See whether the "
                     f"gain plateaus in the table above before pushing higher.\n")
    elif abs(delta) < 1.0:
        lines.append("The effect on convergence is within noise on this sample "
                     f"(Δ < 1 pt). No value materially beats ncrithyst=0; any small "
                     f"non-zero default (e.g. the historical 0.2) is defensible but "
                     f"unsupported by a convergence argument. Iteration counts and "
                     f"the rescue/regression tables above show the practical picture.\n")
    else:
        lines.append("Larger ncrithyst hurt convergence on this sample; prefer a small "
                     "value or 0.\n")

    (OUT_DIR / "findings.md").write_text("\n".join(lines))
    print(f"Wrote {OUT_DIR / 'findings.md'}")


# -------------------------------- main ------------------------------------- #
def main():
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyze-only", action="store_true",
                    help="skip the sweep, analyze existing ncrithyst_raw.csv")
    args = ap.parse_args()

    raw = OUT_DIR / "ncrithyst_raw.csv"
    if args.analyze_only:
        df = pd.read_csv(raw)
    else:
        df = run_sweep()

    stats = analyze(df)
    make_plots(stats)
    write_findings(stats)
    print("Done.")


if __name__ == "__main__":
    main()
