"""Forced-transition AD-vs-FD relative-error vs step-size sweep.
Saves to bench/results/sweep_<label>.npz so before/after can be overlaid.
Usage: python3 bench/grad_sweep.py <label>
"""
import sys
import numpy as np
from GFoil import fwd_run, grad_run, Aerofoil, Acoustics, OperatingConds
from svd301 import doSVD

label = sys.argv[1] if len(sys.argv) > 1 else "after"
xb, yb, Uscaled = doSVD()
nModes = 10
modes  = Uscaled[:, :nModes]
CHORD, SPAN, RE, MA, NCRIT = 0.3, 1.5, 2e6, 0.0, 6.0
OBS = np.array([0.0, 0.0, 3.0]); MODEL='kam'
w0, alf0 = np.zeros(nModes), 3.0
TRANS = [0.1, 0.1]

def build(w, a):
    foil=Aerofoil(xb, yb+modes@w, chord=CHORD, span=SPAN)
    noise=Acoustics(observerXYZ=OBS, model=MODEL, TESampleLoc=0.98)
    op=OperatingConds(alpha=a, Re=RE, Ma=MA, nCrit=NCRIT, transition=np.array(TRANS))
    return foil,noise,op

def ev(w,a):
    foil,noise,op=build(w,a); r=fwd_run(foil,op,noise)
    if not r.converged: return np.nan,np.nan,np.nan
    return r.CL,r.CD,r.OASPL

foil,noise,op=build(w0,alf0); r=fwd_run(foil,op,noise)
g=grad_run(r,foil,op,noise)
adCL=(g.dCL_dy@modes); adCD=(g.dCD_dy@modes); adOA=(g.dOASPL_dy@modes)
adCLa, adCDa, adOAa = g.dCL_dalpha, g.dCD_dalpha, g.dOASPL_dalpha

steps=np.logspace(-3,-8,30)
# rows 0..nModes-1 = modes, row nModes = alpha
errCL=np.zeros((nModes+1, len(steps)))
errCD=np.zeros((nModes+1, len(steps)))
errOA=np.zeros((nModes+1, len(steps)))

def rel(fd, ad): return 100*abs(fd-ad)/max(abs(ad),1e-30)

for s,h in enumerate(steps):
    for m in range(nModes):
        wp=w0.copy(); wp[m]+=h; wn=w0.copy(); wn[m]-=h
        clp,cdp,oap=ev(wp,alf0); cln,cdn,oan=ev(wn,alf0)
        errCL[m,s]=rel((clp-cln)/(2*h), adCL[m])
        errCD[m,s]=rel((cdp-cdn)/(2*h), adCD[m])
        errOA[m,s]=rel((oap-oan)/(2*h), adOA[m])
    clp,cdp,oap=ev(w0,alf0+h); cln,cdn,oan=ev(w0,alf0-h)
    errCL[nModes,s]=rel((clp-cln)/(2*h), adCLa)
    errCD[nModes,s]=rel((cdp-cdn)/(2*h), adCDa)
    errOA[nModes,s]=rel((oap-oan)/(2*h), adOAa)
    print(f"h={h:.2e} done")

np.savez(f"bench/results/sweep_{label}.npz", steps=steps,
         errCL=errCL, errCD=errCD, errOA=errOA,
         adCL=adCL, adCD=adCD, adOA=adOA)
print(f"saved bench/results/sweep_{label}.npz")
