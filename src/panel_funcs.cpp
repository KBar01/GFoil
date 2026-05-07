// Kept separate from panel_funcs.hpp because inviscid_velocity and dvelocity_dgamma
// take the non-template Foil struct (data_structs.h), which differs from the
// template Foil<Real> used in srcAD. Merging would require resolving the
// Isol/Foil struct split — deferred to a later phase.

#include "real_type.h"
#include "data_structs.h"
#include "panel_funcs.hpp"

void inviscid_velocity(const Foil& foil, const Real* gammas, const Real& Vinf,
    const Real& alpha, const Real& CPx, const Real& CPy, Real* velocity)
{
    PanelInfo<Real> info;

    Real a1, b1, a2, b2;
    for (int j = 0; j < Ncoords - 1; ++j) {
        panel_linvortex_velocity<Real>(
            foil.x[colMajorIndex(0,j,2)],
            foil.x[colMajorIndex(1,j,2)],
            foil.x[colMajorIndex(0,j+1,2)],
            foil.x[colMajorIndex(1,j+1,2)],
            CPx, CPy, info, a1, b1, a2, b2);

        velocity[0] += a1*gammas[j] + b1*gammas[j+1];
        velocity[1] += a2*gammas[j] + b2*gammas[j+1];
    }

    panel_constsource_velocity<Real>(
        foil.x[colMajorIndex(0,Ncoords-1,2)],
        foil.x[colMajorIndex(1,Ncoords-1,2)],
        foil.x[colMajorIndex(0,0,2)],
        foil.x[colMajorIndex(1,0,2)],
        CPx, CPy, info, a1, a2);

    Real f1x = a1*(-0.5*foil.te.tcp);
    Real f1y = a2*(-0.5*foil.te.tcp);
    Real f2x = a1*0.5*foil.te.tcp;
    Real f2y = a2*0.5*foil.te.tcp;

    velocity[0] += f1x*gammas[0] + f2x*gammas[Ncoords-1];
    velocity[1] += f1y*gammas[0] + f2y*gammas[Ncoords-1];

    panel_linvortex_velocity<Real>(
        foil.x[colMajorIndex(0,Ncoords-1,2)],
        foil.x[colMajorIndex(1,Ncoords-1,2)],
        foil.x[colMajorIndex(0,0,2)],
        foil.x[colMajorIndex(1,0,2)],
        CPx, CPy, info, a1, b1, a2, b2);

    f1x = (a1+b1)*(0.5*foil.te.tdp);
    f1y = (a2+b2)*(0.5*foil.te.tdp);
    f2x = (a1+b1)*(-0.5*foil.te.tdp);
    f2y = (a2+b2)*(-0.5*foil.te.tdp);

    velocity[0] += f1x*gammas[0] + f2x*gammas[Ncoords-1];
    velocity[1] += f1y*gammas[0] + f2y*gammas[Ncoords-1];

    velocity[0] += Vinf*std::cos(alpha);
    velocity[1] += Vinf*std::sin(alpha);
}


void dvelocity_dgamma(const Foil& foil, const Real& CPx, const Real& CPy, Real* V_G)
{
    PanelInfo<Real> info;

    Real a1, b1, a2, b2;
    for (int j = 0; j < Ncoords - 1; ++j) {
        panel_linvortex_velocity<Real>(
            foil.x[colMajorIndex(0,j,2)],
            foil.x[colMajorIndex(1,j,2)],
            foil.x[colMajorIndex(0,j+1,2)],
            foil.x[colMajorIndex(1,j+1,2)],
            CPx, CPy, info, a1, b1, a2, b2);

        V_G[colMajorIndex(0,j,2)]   += a1;
        V_G[colMajorIndex(1,j,2)]   += a2;
        V_G[colMajorIndex(0,j+1,2)] += b1;
        V_G[colMajorIndex(1,j+1,2)] += b2;
    }

    panel_constsource_velocity<Real>(
        foil.x[colMajorIndex(0,Ncoords-1,2)],
        foil.x[colMajorIndex(1,Ncoords-1,2)],
        foil.x[colMajorIndex(0,0,2)],
        foil.x[colMajorIndex(1,0,2)],
        CPx, CPy, info, a1, a2);

    Real f1x = a1*(-0.5*foil.te.tcp);
    Real f1y = a2*(-0.5*foil.te.tcp);
    Real f2x = a1*0.5*foil.te.tcp;
    Real f2y = a2*0.5*foil.te.tcp;

    V_G[colMajorIndex(0,0,2)]         += f1x;
    V_G[colMajorIndex(1,0,2)]         += f1y;
    V_G[colMajorIndex(0,Ncoords-1,2)] += f2x;
    V_G[colMajorIndex(1,Ncoords-1,2)] += f2y;

    panel_linvortex_velocity<Real>(
        foil.x[colMajorIndex(0,Ncoords-1,2)],
        foil.x[colMajorIndex(1,Ncoords-1,2)],
        foil.x[colMajorIndex(0,0,2)],
        foil.x[colMajorIndex(1,0,2)],
        CPx, CPy, info, a1, b1, a2, b2);

    f1x = (a1+b1)*(0.5*foil.te.tdp);
    f1y = (a2+b2)*(0.5*foil.te.tdp);
    f2x = (a1+b1)*(-0.5*foil.te.tdp);
    f2y = (a2+b2)*(-0.5*foil.te.tdp);

    V_G[colMajorIndex(0,0,2)]         += f1x;
    V_G[colMajorIndex(1,0,2)]         += f1y;
    V_G[colMajorIndex(0,Ncoords-1,2)] += f2x;
    V_G[colMajorIndex(1,Ncoords-1,2)] += f2y;
}
