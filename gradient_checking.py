import numpy as np
from GFoil import fwd_run, grad_run, Aerofoil, Acoustics, OperatingConds
from svd301 import doSVD
from matplotlib import pyplot as plt

# ── case definition ───────────────────────────────────────────────────────────
# Single source of truth for the operating point / geometry config so the AD
# eval and every FD eval are guaranteed identical.
xb, yb, Uscaled = doSVD()
nModes = 10
modes  = Uscaled[:, :nModes]

CHORD = 0.3
SPAN  = 1.5
RE    = 2e6
MA    = 0.0
NCRIT = 6.0
OBS   = np.array([0.0, 0.0, 3.0])   # global freestream-aligned frame, origin at quarter-chord
MODEL = 'kam'

w0   = np.zeros(nModes)
alf0 = 3.0


def _build(weights, alpha):
    y     = yb + modes @ weights
    foil  = Aerofoil(xb, y, chord=CHORD, span=SPAN)
    noise = Acoustics(observerXYZ=OBS, model=MODEL,TESampleLoc=0.98)
    op    = OperatingConds(alpha=alpha, Re=RE, Ma=MA, nCrit=NCRIT,transition=np.array([0.1,0.1]))
    return foil, noise, op


def evalPoint(weights, alpha):
    foil, noise, op = _build(weights, alpha)
    result = fwd_run(foil, op, noise)
    if not result.converged:
        # Returning NaNs makes a non-converged perturbation show up as a NaN in
        # the FD curve rather than silently injecting a stale/garbage value.
        raise RuntimeError(
            f"fwd_run did not converge (alpha={alpha}, "
            f"failure_mode={result.failure_mode})"
        )
    return result.CL, result.CD, result.OASPL


def _safe_eval(weights, alpha):
    try:
        return evalPoint(weights, alpha)
    except RuntimeError as e:
        print(f"  WARN: {e}")
        return np.nan, np.nan, np.nan


def singleFD(weights, alpha, step):
    dCLdw    = np.zeros_like(weights)
    dCDdw    = np.zeros_like(weights)
    dOASPLdw = np.zeros_like(weights)

    for i in range(weights.shape[0]):
        wp = weights.copy(); wp[i] += step
        wn = weights.copy(); wn[i] -= step

        CLp, CDp, OASPLp = _safe_eval(wp, alpha)
        CLn, CDn, OASPLn = _safe_eval(wn, alpha)

        dCLdw[i]    = (CLp - CLn)       / (2 * step)
        dCDdw[i]    = (CDp - CDn)       / (2 * step)
        dOASPLdw[i] = (OASPLp - OASPLn) / (2 * step)

    CLp, CDp, OASPLp = _safe_eval(weights, alpha + step)
    CLn, CDn, OASPLn = _safe_eval(weights, alpha - step)
    dCLda    = (CLp - CLn)       / (2 * step)
    dCDda    = (CDp - CDn)       / (2 * step)
    dOASPLda = (OASPLp - OASPLn) / (2 * step)

    return dCLdw, dCLda, dCDdw, dCDda, dOASPLdw, dOASPLda


# ── AD gradients ──────────────────────────────────────────────────────────────
foil, noise, op = _build(w0, alf0)

# Guard against silent surface-ordering flips in Aerofoil.__post_init__. If GFoil
# reverses the coordinate arrays to normalise orientation, the row ordering of
# `modes` no longer matches dCx_dy and the chain rule below would be wrong.
y0 = yb + modes @ w0
if not np.allclose(foil.ycoords, y0):
    raise RuntimeError(
        "Aerofoil normalised/flipped the coordinate ordering; the SVD `modes` "
        "basis is no longer aligned with dCx_dy. Re-order `modes` (and xb/yb) "
        "to match foil.xcoords/foil.ycoords before contracting."
    )

result = fwd_run(foil, op, noise)
if not result.converged:
    raise RuntimeError(f"Baseline fwd_run did not converge: {result.failure_mode}")
print("done init eval ---------------")
print(result.summary())

grads = grad_run(result, foil, op, noise)
if not grads.converged:
    raise RuntimeError("grad_run reported non-converged forward state.")

dCLdy    = grads.dCL_dy
dCDdy    = grads.dCD_dy
dNoisedy = grads.dOASPL_dy

dCLdalpha    = grads.dCL_dalpha
dCDdalpha    = grads.dCD_dalpha
dOASPLdalpha = grads.dOASPL_dalpha

# project coordinate sensitivities onto the SVD design modes
dCLdweights    = dCLdy    @ modes
dCDdweights    = dCDdy    @ modes
dOASPLdweights = dNoisedy @ modes

np.savetxt('AD_OASPL_weights.txt', np.asarray(dOASPLdweights), delimiter=',')
np.savetxt('AD_CD_weights.txt',    np.asarray(dCDdweights),    delimiter=',')
np.savetxt('AD_CL_weights.txt',    np.asarray(dCLdweights),    delimiter=',')

np.savetxt('AD_OASPL_alpha.txt', np.atleast_1d(dOASPLdalpha))
np.savetxt('AD_CD_alpha.txt',    np.atleast_1d(dCDdalpha))
np.savetxt('AD_CL_alpha.txt',    np.atleast_1d(dCLdalpha))

# ── finite difference sweep ───────────────────────────────────────────────────
steps = np.logspace(-3, -8, 50)

CLgrads    = np.zeros([nModes + 1, steps.shape[0]])
CDgrads    = np.zeros([nModes + 1, steps.shape[0]])
OASPLgrads = np.zeros([nModes + 1, steps.shape[0]])

for s in range(steps.shape[0]):
    print(f"step {s} out of {steps.shape[0]}")
    h = steps[s]
    dCLdw, dCLda, dCDdw, dCDda, dOASPLdw, dOASPLda = singleFD(w0, alf0, h)
    CLgrads[:-1, s]    = dCLdw
    CLgrads[-1,  s]    = dCLda
    CDgrads[:-1, s]    = dCDdw
    CDgrads[-1,  s]    = dCDda
    OASPLgrads[:-1, s] = dOASPLdw
    OASPLgrads[-1,  s] = dOASPLda

np.savetxt('FD_OASPL_weights.txt', np.asarray(OASPLgrads), delimiter=',')
np.savetxt('FD_CD_weights.txt',    np.asarray(CDgrads),    delimiter=',')
np.savetxt('FD_CL_weights.txt',    np.asarray(CLgrads),    delimiter=',')

# ── plots ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "mathtext.fontset": "stix",
    "font.family":      "STIXGeneral",
    "font.size":        14,
    "axes.labelsize":   14,
    "axes.titlesize":   14,
    "legend.fontsize":  14,
    "xtick.labelsize":  14,
    "ytick.labelsize":  14,
})


def _rel_err(fd_row, ad_val):
    # percentage relative error, guarding a zero/near-zero AD reference
    denom = np.abs(ad_val)
    if denom < 1e-30:
        return np.full_like(fd_row, np.nan)
    return 100.0 * np.abs(fd_row - ad_val) / denom


def make_plot(fd_grads, ad_weights, ad_alpha, ylabel, fname):
    plt.figure()
    for i in range(nModes):
        plt.loglog(steps, _rel_err(fd_grads[i, :], ad_weights[i]), label=f"mode {i+1}")
    plt.loglog(steps, _rel_err(fd_grads[-1, :], ad_alpha), label="alpha")
    plt.xlabel('step size')
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(fname, dpi=300)
    plt.close()


make_plot(CLgrads,    dCLdweights,    dCLdalpha,    'CL error, %',    "CLerror.png")
make_plot(CDgrads,    dCDdweights,    dCDdalpha,    'CD error, %',    "CDerror.png")
make_plot(OASPLgrads, dOASPLdweights, dOASPLdalpha, 'OASPL error, %', "OASPLerror.png")