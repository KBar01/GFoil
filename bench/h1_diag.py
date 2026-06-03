"""H1 diagnostic: AD-vs-FD agreement, FORCED vs FREE transition.

Hypothesis H1 predicts the forced-transition arc-length station xift is detached
from the CoDi tape (computed via .getValue()), so d(xift)/dy is missing from the
adjoint. If true: AD-FD agreement should be good in FREE transition but broken in
FORCED transition, and the error should scale with how much each SVD mode displaces
the surface near x/c = 0.1.
"""
import sys
import numpy as np
from GFoil import fwd_run, grad_run, Aerofoil, Acoustics, OperatingConds
from svd301 import doSVD

xb, yb, Uscaled = doSVD()
nModes = 10
modes  = Uscaled[:, :nModes]

CHORD, SPAN, RE, MA, NCRIT = 0.3, 1.5, 2e6, 0.0, 6.0
OBS   = np.array([0.0, 0.0, 3.0])
MODEL = 'kam'
w0, alf0 = np.zeros(nModes), 3.0


def build(weights, alpha, trans):
    y    = yb + modes @ weights
    foil = Aerofoil(xb, y, chord=CHORD, span=SPAN)
    noise= Acoustics(observerXYZ=OBS, model=MODEL, TESampleLoc=0.98)
    op   = OperatingConds(alpha=alpha, Re=RE, Ma=MA, nCrit=NCRIT,
                          transition=np.array(trans))
    return foil, noise, op


def evalP(weights, alpha, trans):
    foil, noise, op = build(weights, alpha, trans)
    r = fwd_run(foil, op, noise)
    if not r.converged:
        raise RuntimeError(f"no converge a={alpha} mode-perturb fm={r.failure_mode}")
    return r.CL, r.CD, r.OASPL


def ad_grads(trans):
    foil, noise, op = build(w0, alf0, trans)
    r = fwd_run(foil, op, noise)
    assert r.converged
    g = grad_run(r, foil, op, noise)
    return (g.dCL_dy @ modes, g.dCD_dy @ modes, g.dOASPL_dy @ modes,
            g.dCL_dalpha, g.dCD_dalpha, g.dOASPL_dalpha, r.OASPL)


def fd_col(weights, alpha, trans, h):
    dCL, dCD, dOA = np.zeros(nModes), np.zeros(nModes), np.zeros(nModes)
    for i in range(nModes):
        wp = weights.copy(); wp[i] += h
        wn = weights.copy(); wn[i] -= h
        clp, cdp, oap = evalP(wp, alpha, trans)
        cln, cdn, oan = evalP(wn, alpha, trans)
        dCL[i] = (clp-cln)/(2*h); dCD[i] = (cdp-cdn)/(2*h); dOA[i] = (oap-oan)/(2*h)
    clp, cdp, oap = evalP(weights, alpha+h, trans)
    cln, cdn, oan = evalP(weights, alpha-h, trans)
    da = ((clp-cln)/(2*h), (cdp-cdn)/(2*h), (oap-oan)/(2*h))
    return dCL, dCD, dOA, da


def relerr(fd, ad):
    return 100.0*np.abs(fd-ad)/np.maximum(np.abs(ad), 1e-30)


def run(label, trans, h):
    print(f"\n===== {label}  (transition={trans}, h={h}) =====")
    adCL, adCD, adOA, adCLa, adCDa, adOAa, oaspl = ad_grads(trans)
    print(f"baseline OASPL = {oaspl:.6f} dB")
    fdCL, fdCD, fdOA, (fdCLa, fdCDa, fdOAa) = fd_col(w0, alf0, trans, h)
    # surface displacement near x/c=0.1 per mode (proxy for H1 scaling)
    xc = xb / xb.max()
    near = np.abs(xc - 0.1) < 0.03
    mode_disp = np.array([np.max(np.abs(modes[near, i])) for i in range(nModes)])
    print(f"{'mode':>5} {'|disp@0.1|':>11} {'OASPL_AD':>11} {'OASPL_FD':>11} "
          f"{'OASPL%':>9} {'CL%':>9} {'CD%':>9}")
    for i in range(nModes):
        print(f"{i+1:>5} {mode_disp[i]:>11.3e} {adOA[i]:>11.4e} {fdOA[i]:>11.4e} "
              f"{relerr(fdOA[i],adOA[i]):>9.3f} {relerr(fdCL[i],adCL[i]):>9.3f} "
              f"{relerr(fdCD[i],adCD[i]):>9.3f}")
    print(f"alpha: OASPL AD={adOAa:.5e} FD={fdOAa:.5e} err={relerr(fdOAa,adOAa):.3f}%  "
          f"| CL AD={adCLa:.5e} FD={fdCLa:.5e} err={relerr(fdCLa,adCLa):.3f}%  "
          f"| CD AD={adCDa:.5e} FD={fdCDa:.5e} err={relerr(fdCDa,adCDa):.3f}%")
    return relerr(fdOA, adOA), mode_disp


if __name__ == "__main__":
    h = 1e-5
    free_err, _      = run("FREE transition",   [1.0, 1.0], h)
    forced_err, disp = run("FORCED transition", [0.1, 0.1], h)
    print("\n----- H1 summary (OASPL mode-error, forced vs free) -----")
    print(f"{'mode':>5} {'|disp@0.1|':>11} {'free%':>10} {'forced%':>10}")
    for i in range(nModes):
        print(f"{i+1:>5} {disp[i]:>11.3e} {free_err[i]:>10.3f} {forced_err[i]:>10.3f}")
