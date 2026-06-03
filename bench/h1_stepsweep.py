"""Step-size sweep, FORCED transition. Distinguish:
  H1 signature  -> flat error FLOOR independent of h (constant missing term)
  H2 signature  -> error CLIFF: high for large h, sudden drop when h small enough
                   to stop straddling a discontinuity.
"""
import numpy as np
from GFoil import fwd_run, grad_run, Aerofoil, Acoustics, OperatingConds
from svd301 import doSVD

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
adOA=(g.dOASPL_dy@modes); adCD=(g.dCD_dy@modes); adCL=(g.dCL_dy@modes)

probe=[0,1,3,6,8]   # modes 1,2,4,7,9
steps=np.logspace(-3,-8,20)
print("OASPL relative error (%) vs step, forced transition")
print("h        " + "".join(f" mode{m+1:>2}" for m in probe))
for h in steps:
    row=[]
    for m in probe:
        wp=w0.copy(); wp[m]+=h; wn=w0.copy(); wn[m]-=h
        _,_,op_=ev(wp,alf0); _,_,on_=ev(wn,alf0)
        fd=(op_-on_)/(2*h)
        row.append(100*abs(fd-adOA[m])/max(abs(adOA[m]),1e-30))
    print(f"{h:.2e} " + "".join(f"{v:8.2f}" for v in row))
