#pragma once
#include "real_type.hpp"
#include "data_structs.hpp"
#include "vector_ops.hpp"
#include "main_func.hpp"   // identify_surfaces
#include "calc_ue_m.hpp"   // rebuild_ue_m

template<typename Real>
void stagpoint_move(Isolv<Real>& isol, Glob<Real>& glob, const Foil<Real>& foil,
    const Wake<Real>& wake, Vsol<Real>& vsol)
{
    int* I = isol.stagIndex;
    bool newpanel = true;
    const int nc = foil.nc, nw = wake.nw;

    if (glob.U[colMajorIndex(3,I[1],4)] < 0) {
        int j = I[1];
        for (; j < nc; ++j) {
            if (glob.U[colMajorIndex(3,j,4)] > 0) break;
        }
        int I1 = j;
        for (j = I[1]; j < I1; ++j) glob.U[colMajorIndex(3,j,4)] *= -1.0;
        I[0] = I1-1; I[1] = I1;
    } else if (glob.U[colMajorIndex(3,I[0],4)] < 0) {
        int j = I[0];
        for (; j >= 0; --j) {
            if (glob.U[colMajorIndex(3,j,4)] > 0) break;
        }
        assert(j > 0 && "no stagnation point");
        int I0 = j;
        for (j = I0+1; j <= I[0]; ++j) glob.U[colMajorIndex(3,j,4)] *= -1.0;
        I[0] = I0; I[1] = I0+1;
    } else {
        newpanel = false;
    }

    Real u0 = glob.U[colMajorIndex(3,I[0],4)];
    Real u1 = glob.U[colMajorIndex(3,I[1],4)];
    Real den = u0+u1;
    Real w1 = u1/den, w2 = u0/den;

    isol.stagArcLocation = w1*foil.s[I[0]] + w2*foil.s[I[1]];
    for (int d=0;d<2;++d)
        isol.stagXLocation[d] = w1*foil.x[colMajorIndex(d,I[0],2)] + w2*foil.x[colMajorIndex(d,I[1],2)];

    Real ds = foil.s[I[1]]-foil.s[I[0]];
    isol.sstag_ue[0] =  u1*ds/(den*den);
    isol.sstag_ue[1] = -u0*ds/(den*den);

    for (int i=0; i<nc; ++i)
        isol.distFromStag[i] = std::abs(foil.s[i] - isol.stagArcLocation);
    for (int i=0; i<nw; ++i)
        isol.distFromStag[nc+i] = wake.s[i] - isol.stagArcLocation;

    if (newpanel) {
        for (int i=0; i<=I[0]; ++i) isol.edgeVelSign[i] = -1;
        for (int i=I[0]+1; i<nc; ++i) isol.edgeVelSign[i] = 1;
        identify_surfaces<Real>(isol, vsol);
        rebuild_ue_m<Real>(foil, wake, isol, vsol, true);
    }
}
