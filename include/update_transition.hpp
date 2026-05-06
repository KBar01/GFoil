#pragma once
#include "real_type.hpp"
#include "data_structs.hpp"
#include "get_funcs.hpp"
#include "residuals.hpp"

template<typename Real>
int march_amplification(Glob<Real>& glob, Vsol<Real>& vsol, const Isolv<Real>& isol,
    int si, const Param<Real>& param, Trans<Real>& tdata, const bool force)
{
    const std::vector<int>& Is = vsol.Is[si];
    int N = Is.size();
    glob.U[colMajorIndex(2,Is[0],4)] = 0.0;

    int i = 1;
    while (i < N) {
        int i1 = Is[i-1], i2 = Is[i];
        Real U1[4], U2[4];
        for (int j=0;j<4;++j) {
            U1[j] = glob.U[colMajorIndex(j,i1,4)];
            U2[j] = glob.U[colMajorIndex(j,i2,4)];
        }
        if (vsol.turb[i2]) U2[2] = U1[2]*1.01;

        Real dx = isol.distFromStag[i2]-isol.distFromStag[i1];
        constexpr int nNewton = 20;
        const Real one=1;

        for (int iNewton=0; iNewton<nNewton; ++iNewton) {
            Real damp1_U[4], damp2_U[4], damp_U[8]={0};
            Real damp1 = get_damp(U1[0],U1[1],U1[2],U1[3],param,damp1_U);
            Real damp2 = get_damp(U2[0],U2[1],U2[2],U2[3],param,damp2_U);
            Real damp  = upwind_half(damp1,damp1_U,damp2,damp2_U,damp_U);

            Real Ramp = U2[2]-U1[2]-damp*dx;
            if (std::abs(Ramp)<1e-12) break;

            Real Ramp_U[8]={0};
            Ramp_U[2]=-1; Ramp_U[6]=1;
            for (int j=0;j<8;++j) Ramp_U[j] -= damp_U[j]*dx;

            Real dU = -Ramp/Ramp_U[6];
            Real dmax = 0.5*(1.01-Real(iNewton)/nNewton);
            Real om = (std::abs(dU)>dmax) ? dmax/std::abs(dU) : one;
            U2[2] += om*dU;
        }

        if (U2[2] > param.ncrit) { tdata.isForced[si]=0; break; }
        else if (force && i1==tdata.transNode[si]) { tdata.isForced[si]=1; break; }
        else { glob.U[colMajorIndex(2,i2,4)] = U2[2]; }
        ++i;
    }
    return i-1;
}

template<typename Real>
void update_transition(Glob<Real>& glob, Vsol<Real>& vsol, const Isolv<Real>& isol,
    Param<Real>& param, Trans<Real>& tdata, const bool force)
{
    for (int si=0;si<2;++si) {
        const std::vector<int>& Is = vsol.Is[si];
        int nSurfPoints = Is.size();

        int ilam0 = nSurfPoints-1;
        for (int i=0;i<nSurfPoints;++i) {
            if (vsol.turb[Is[i]]) { ilam0=i-1; break; }
        }

        std::vector<Real> sa(glob.nc);
        for (int s=0;s<glob.nc;++s) sa[s] = glob.U[colMajorIndex(2,s,4)];

        int ilam = march_amplification<Real>(glob,vsol,isol,si,param,tdata,force);

        if (ilam==ilam0) {
            for (int s=0;s<glob.nc;++s) glob.U[colMajorIndex(2,s,4)]=sa[s];
            continue;
        }

        if (ilam < ilam0) {
            Real cttr_U[4];
            Real sa0 = get_cttr(glob.U[colMajorIndex(0,Is[ilam+1],4)],
                glob.U[colMajorIndex(1,Is[ilam+1],4)],
                glob.U[colMajorIndex(2,Is[ilam+1],4)],
                glob.U[colMajorIndex(3,Is[ilam+1],4)],
                true,param,cttr_U);
            Real sa1 = (ilam0<nSurfPoints-1) ? glob.U[colMajorIndex(2,Is[ilam0+1],4)] : sa0;

            const Real zero=0, one=1;
            Real xi_start = isol.distFromStag[Is[ilam+1]];
            Real xi_end   = isol.distFromStag[Is[std::min(ilam0+1,nSurfPoints-1)]];
            Real dx = xi_end-xi_start;

            for (int i=ilam+1; i<=ilam0; ++i) {
                Real f = (dx==0||i==ilam+1) ? zero : (isol.distFromStag[Is[i]]-xi_start)/dx;
                if ((ilam+1)==ilam0) f=one;
                glob.U[colMajorIndex(2,Is[i],4)] = sa0+f*(sa1-sa0);
                vsol.turb[Is[i]] = 1;
            }
        } else {
            for (int i=ilam0; i<=ilam; ++i) vsol.turb[Is[i]]=0;
        }
    }
}
