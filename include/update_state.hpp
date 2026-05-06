#pragma once
#include "real_type.hpp"
#include "data_structs.hpp"
#include "get_funcs.hpp"

template<typename Real>
void update_state(const Oper<Real>& oper, const Param<Real>& param, Glob<Real>& glob, Vsol<Real>& vsol)
{
    const int Nsys = glob.nc + glob.nw;

    std::vector<int> It;
    Real ctmax = -1e20;
    for (int i=0;i<Nsys;++i) {
        if (vsol.turb[i]) {
            It.push_back(i);
            Real ct = glob.U[colMajorIndex(2,i,4)];
            if (ct > ctmax) ctmax = ct;
        }
    }

    Real omega = 1.0;
    const Real one = 1;

    for (int k=0;k<2;++k) {
        Real fmin = 1e20;
        for (int i=0;i<Nsys;++i) {
            Real Uk  = glob.U[colMajorIndex(k,i,4)];
            Real dUk = glob.dU[colMajorIndex(k,i,4)];
            if (Uk != 0.0) { Real r = dUk/Uk; if (r<fmin) fmin=r; }
        }
        Real om = (fmin < -0.5) ? std::abs(0.5/fmin) : one;
        if (om<omega) omega=om;
    }

    for (int i=0;i<Nsys;++i) {
        Real Uk  = glob.U[colMajorIndex(2,i,4)];
        Real dUk = glob.dU[colMajorIndex(2,i,4)];
        if (!vsol.turb[i] && Uk<0.2) continue;
        if ( vsol.turb[i] && Uk<0.1*ctmax) continue;
        if (Uk==0.0 || dUk==0.0) continue;
        if ((Uk+dUk)<0.0) { Real om=0.8*std::abs(Uk/dUk); if (om<omega) omega=om; }
    }

    Real dumax=0;
    for (int i=0;i<Nsys;++i) {
        if (!vsol.turb[i]) { Real d=std::abs(glob.dU[colMajorIndex(2,i,4)]); if (d>dumax) dumax=d; }
    }
    { Real om=(dumax>0)?std::abs(2.0/dumax):one; if (om<omega) omega=om; }

    dumax=0;
    for (int i=0;i<Nsys;++i) {
        if (vsol.turb[i]) { Real d=std::abs(glob.dU[colMajorIndex(2,i,4)]); if (d>dumax) dumax=d; }
    }
    { Real om=(dumax>0)?std::abs(0.05/dumax):one; if (om<omega) omega=om; }

    dumax=0;
    for (int i=0;i<Nsys;++i) {
        Real d=std::abs(glob.dU[colMajorIndex(3,i,4)]/oper.Vinf);
        if (d>dumax) dumax=d;
    }
    { Real om=(dumax>0)?0.2/dumax:one; if (om<omega) omega=om; }

    for (int i=0;i<4*Nsys;++i) glob.U[i] += omega*glob.dU[i];

    for (int si=0;si<3;++si) {
        Real Hkmin = (si==2) ? 1.00005 : 1.02;
        const std::vector<int>& Is = vsol.Is[si];
        for (int j : Is) {
            Real Uj[4] = {glob.U[colMajorIndex(0,j,4)],glob.U[colMajorIndex(1,j,4)],
                          glob.U[colMajorIndex(2,j,4)],glob.U[colMajorIndex(3,j,4)]};
            Real Hk_U[4];
            Real Hk = get_Hk(Uj[0],Uj[1],Uj[3],param,Hk_U);
            if (Hk < Hkmin)
                glob.U[colMajorIndex(1,j,4)] += 2.0*(Hkmin-Hk)*glob.U[colMajorIndex(1,j,4)];
        }
    }

    for (int i : It) {
        if (glob.U[colMajorIndex(2,i,4)] < 0.0)
            glob.U[colMajorIndex(2,i,4)] = 0.1*ctmax;
    }
}
