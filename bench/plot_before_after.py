"""Before/after AD-vs-FD log-log error plots, forced transition.
Reads bench/results/sweep_{before,after}.npz, writes 3 PNGs to bench/results/."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = np.load("bench/results/sweep_before.npz")
A = np.load("bench/results/sweep_after.npz")
steps = B["steps"]
nModes = B["errOA"].shape[0] - 1

def plot(key, ylabel, fname):
    fig, ax = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for col, (D, ttl) in enumerate([(B, "BEFORE (xift detached)"),
                                     (A, "AFTER (xift taped)")]):
        E = D[key]
        for i in range(nModes):
            ax[col].loglog(steps, E[i], label=f"mode {i+1}", lw=1)
        ax[col].loglog(steps, E[nModes], "k--", lw=2, label="alpha")
        ax[col].set_title(ttl)
        ax[col].set_xlabel("FD step size h")
        ax[col].grid(True, which="both", alpha=0.3)
    ax[0].set_ylabel(ylabel)
    ax[1].legend(fontsize=8, ncol=2, loc="upper right")
    fig.suptitle(f"{ylabel} — forced transition (x/c=0.1), NACA0012 SVD foil, "
                 f"Re=2e6, alpha=3deg")
    fig.tight_layout()
    fig.savefig(f"bench/results/{fname}", dpi=150)
    plt.close(fig)
    print(f"wrote bench/results/{fname}")

plot("errCL", "CL rel. error, %",    "GRAD_VERIFY_CL.png")
plot("errCD", "CD rel. error, %",    "GRAD_VERIFY_CD.png")
plot("errOA", "OASPL rel. error, %", "GRAD_VERIFY_OASPL.png")

# Print a compact before/after worst-error summary (best step per row).
print("\nWorst-case (min-over-h) rel error %, forced transition:")
print(f"{'row':>6} | {'CL before':>9} {'CL after':>9} | {'CD before':>9} {'CD after':>9} | {'OA before':>9} {'OA after':>9}")
labels = [f"mode{i+1}" for i in range(nModes)] + ["alpha"]
for i in range(nModes + 1):
    print(f"{labels[i]:>6} | {B['errCL'][i].min():9.3f} {A['errCL'][i].min():9.3f} | "
          f"{B['errCD'][i].min():9.3f} {A['errCD'][i].min():9.3f} | "
          f"{B['errOA'][i].min():9.3f} {A['errOA'][i].min():9.3f}")
