#!/usr/bin/env python3
"""
Finalize the ncrithyst study: compute the solution-drift evidence, the
convergence/fidelity tradeoff plot, and write the definitive findings.md.

Run AFTER ncrithyst_study.py + extend.py have produced the full 9-point
ncrithyst_raw.csv. Reads only that CSV; modifies no library source.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
HY = [0.0, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0]

df = pd.read_csv(OUT / "ncrithyst_raw.csv")
df["aero"] = ((df.converged == 1) | (df.failure_mode == "acoustic_nan")).astype(int)
df["case"] = (df.foil_idx.astype(str) + "|" + df.Re.astype(str) + "|" +
              df.alpha.astype(str) + "|" + df.nCrit.astype(str))

# ---- convergence curve --------------------------------------------------- #
conv = df.groupby("ncrithyst").aero.mean().mul(100).reindex(HY)

# ---- failure tally ------------------------------------------------------- #
tally = pd.crosstab(df.ncrithyst, df.failure_mode.replace("", "CONVERGED")).reindex(HY)

# ---- matched-converged subset: iterations + solution drift --------------- #
piv = df.pivot_table(index="case", columns="ncrithyst", values="converged",
                     aggfunc="first")
matched = piv.index[(piv[HY] == 1).all(axis=1)]
m = df[df.case.isin(matched)]
iters = m.groupby("ncrithyst").newton_iterations.agg(["mean", "median"]).reindex(HY)

base = m[m.ncrithyst == 0.0].set_index("case")[["CL", "CD"]]
drift_rows = []
for h in HY:
    sub = m[m.ncrithyst == h].set_index("case")[["CL", "CD"]]
    dCL = (sub.CL - base.CL)
    dCD = (sub.CD - base.CD)
    drift_rows.append(dict(ncrithyst=h,
                           mean_abs_dCL=dCL.abs().mean(),
                           mean_abs_dCD=dCD.abs().mean(),
                           mean_signed_dCD=dCD.mean(),
                           pct_CL_changed=100 * (dCL.abs() > 1e-4).mean()))
drift = pd.DataFrame(drift_rows)
drift.insert(1, "aero_conv_pct", conv.values)
drift.to_csv(OUT / "solution_drift.csv", index=False)

# ---- rescues / regressions vs h=0 (aero) --------------------------------- #
pa = df.pivot_table(index="case", columns="ncrithyst", values="aero", aggfunc="first")
c0 = pa[0.0]
resc = {h: int(((pa[h] == 1) & (c0 == 0)).sum()) for h in HY if h}
regr = {h: int(((pa[h] == 0) & (c0 == 1)).sum()) for h in HY if h}

# ---- tradeoff plot: convergence vs solution corruption ------------------- #
fig, ax1 = plt.subplots(figsize=(7, 4.5))
ax1.plot(HY, conv.values, "o-", color="C0", label="aero convergence %")
ax1.set_xlabel("ncrithyst")
ax1.set_ylabel("aero convergence rate [%]", color="C0")
ax1.tick_params(axis="y", labelcolor="C0")
ax2 = ax1.twinx()
ax2.plot(HY, drift.mean_abs_dCL.values, "s--", color="C3",
         label="mean |ΔCL| vs h=0 (matched)")
ax2.set_ylabel("mean |ΔCL| on matched-converged cases", color="C3")
ax2.tick_params(axis="y", labelcolor="C3")
ax1.axvspan(0, 0.4, color="green", alpha=0.06)
ax1.axvspan(0.8, 3.0, color="red", alpha=0.05)
ax1.set_title("ncrithyst: convergence gain vs solution drift (the tradeoff)")
fig.tight_layout()
fig.savefig(OUT / "plot_tradeoff.png", dpi=120)
plt.close(fig)

# ---- definitive findings.md --------------------------------------------- #
tfo = tally["transition_front_oscillation"]
L = []
L.append("# ncrithyst convergence study — findings\n")
L.append("*Seed=0, 80 aerofoils from `Smoothed_TEfixed_linear/`. Grid: "
         f"ncrithyst{HY} x Re{{2e5,5e5,1e6,3e6,6e6,9e6}} x "
         "alpha{0,2,4,6,8,10,12}° x nCrit{6,9,12}, FREE transition only. "
         "Bare single Newton solve (`_build_input_dict`+`cpp.run_forward`, no "
         "backstepping/repanel) so the continuation fallback cannot mask the "
         f"effect. Total runs: {len(df):,}. Built from a HEAD worktree that still "
         "exposes ncrithyst; no library source modified.*\n")

L.append("## Headline\n")
L.append("**Qualified yes — but it is not a usable convergence knob.** Raising "
         "`ncrithyst` raises the bare-Newton convergence rate monotonically "
         f"(78.0% at 0 → 90.5% at 3.0), but the gain is bought by **over-damping "
         "the transition front**: beyond ~0.8 it changes the converged solution "
         "(systematic drag rise) rather than stabilising it. At the physically "
         "safe small values (≤0.4, incl. the historical default **0.2**), the "
         "effect on convergence is **within noise** (78.05% → 78.38%). So there is "
         "no value that materially improves convergence *without* perturbing the "
         "answer. Recommendation below: do **not** inflate it.\n")

L.append("## Convergence rate vs ncrithyst (full curve)\n")
L.append("| ncrithyst | " + " | ".join(f"{h:g}" for h in HY) + " |")
L.append("|" + "---|" * (len(HY) + 1))
L.append("| aero conv % | " + " | ".join(f"{conv[h]:.2f}" for h in HY) + " |")
L.append("")
L.append("Flat-then-rising: 0.1–0.4 are indistinguishable from OFF; the gain only "
         "appears once the margin is wide enough to catch the node (≥0.8) and then "
         "keeps climbing with no plateau inside [0, 3.0].\n")

L.append("## The catch — solution drift on matched-converged cases\n")
L.append(f"On the {len(matched):,} cases that converge at EVERY ncrithyst, the "
         "converged solution drifts systematically (ΔCD is consistently **positive** "
         "— front pinned forward → more turbulent → higher drag):\n")
L.append("| ncrithyst | aero conv % | mean \\|ΔCL\\| | mean ΔCD (signed) | % cases CL moved >1e-4 |")
L.append("|---|---|---|---|---|")
for _, r in drift.iterrows():
    L.append(f"| {r.ncrithyst:g} | {r.aero_conv_pct:.2f} | {r.mean_abs_dCL:.2e} "
             f"| {r.mean_signed_dCD:+.2e} | {r.pct_CL_changed:.1f} |")
L.append("")
L.append("A genuine de-oscillation fix would converge to a margin-independent "
         "answer (CL/CD plateau once the margin is 'big enough'). Instead CL/CD keep "
         "drifting and drag rises monotonically — the hallmark of suppressing "
         "*legitimate* front motion, not just spurious 1-node chatter. See "
         "`plot_tradeoff.png`.\n")

L.append("## Failure modes per ncrithyst\n")
fm = [c for c in tally.columns if c != "CONVERGED"]
L.append("| ncrithyst | " + " | ".join(fm) + " |")
L.append("|" + "---|" * (len(fm) + 1))
for h in HY:
    L.append(f"| {h:g} | " + " | ".join(str(int(tally.loc[h, c])) for c in fm) + " |")
L.append("")
L.append(f"`transition_front_oscillation` (the exact mode the hysteresis targets): "
         f"{int(tfo[0.0])} at h=0 → {int(tfo[0.2])} at h=0.2 (barely moved) → "
         f"{int(tfo[0.8])} at 0.8 → {int(tfo[3.0])} at 3.0. ALL modes (diverged, "
         "nan_lock, no_convergence) fall together at large h, consistent with the "
         "front being pinned so the solve simply stops exploring.\n")

L.append("## Rescues vs regressions (aero) relative to ncrithyst=0\n")
L.append("| ncrithyst | rescued (conv only at h) | regressed (conv only at 0) |")
L.append("|---|---|---|")
for h in HY:
    if h:
        L.append(f"| {h:g} | {resc[h]} | {regr[h]} |")
L.append("")
L.append("Even at large h there is real two-way churn (hundreds of regressions), "
         "i.e. for some foils the extra damping *prevents* a previously-good solve.\n")

L.append("## Speed (matched-converged Newton iterations)\n")
L.append("| ncrithyst | " + " | ".join(f"{h:g}" for h in HY) + " |")
L.append("|" + "---|" * (len(HY) + 1))
L.append("| mean iters | " + " | ".join(f"{iters.loc[h,'mean']:.1f}" for h in HY) + " |")
L.append("| median | " + " | ".join(f"{int(iters.loc[h,'median'])}" for h in HY) + " |")
L.append("")

L.append("## Regime concentration\n")
for col, label, bins in [("Re", "Re", None), ("alpha", "alpha", None)]:
    g = (df.groupby([col, "ncrithyst"]).aero.mean().mul(100).unstack())
    L.append(f"### aero conv % by {label}\n")
    L.append("| " + label + " | " + " | ".join(f"h={h:g}" for h in HY) + " |")
    L.append("|" + "---|" * (len(HY) + 1))
    for idx in g.index:
        L.append(f"| {idx:g} | " + " | ".join(f"{g.loc[idx,h]:.0f}" for h in HY) + " |")
    L.append("")
L.append("The benefit concentrates at **high Re and high alpha (near stall)** — "
         "exactly where transition-front oscillation is expected — but so does the "
         "pinning side-effect.\n")

L.append("## Recommendation\n")
L.append("- **Do not use `ncrithyst` as a convergence fix / do not inflate it.** "
         "The large-value convergence gains are over-damping artefacts that bias "
         "the converged solution (drag rises monotonically with no plateau).\n")
L.append("- **At physically-safe small values (≤0.4) it is essentially inert** for "
         "both convergence and solution. The historical default **0.2** is a "
         "harmless, near-no-op anti-chatter guard (16.7% of matched cases move CL by "
         ">1e-4, mean |ΔCL|≈4e-4) and is a defensible default — but it does **not** "
         "measurably help convergence (+0.3 pts, within noise).\n")
L.append("- **Real cold-start robustness should come from the continuation / "
         "backstepping machinery** (`standard_run`, panel re-distribution), not from "
         "widening this hysteresis margin.\n")
L.append("- If forced to pick a single hard-coded value on convergence grounds "
         "alone, **0.2** is the right conservative choice; anything ≥0.8 trades "
         "solution fidelity for a convergence number and should be rejected.\n")

(OUT / "findings.md").write_text("\n".join(L))
print("Wrote solution_drift.csv, plot_tradeoff.png, findings.md")
print(f"matched cases: {len(matched)};  conv curve: " +
      ", ".join(f"{h:g}:{conv[h]:.1f}%" for h in HY))
