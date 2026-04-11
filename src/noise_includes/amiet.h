#ifndef AMIET_H
#define AMIET_H

#include "real_type.h"


void calc_WPS_Goody(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real edgeVel,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real Uinf,
                    Real (&phiqq)[Nsound]
);


void calc_WPS_Kamruzzaman(Real theta,
                    Real deltaS,
                    Real tauWall,
                    Real Ue,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real (&phiqq)[Nsound]);

void calc_WPS_Rozenburg(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real Ue,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real (&phiqq)[Nsound]
);

void calc_WPS_Lee(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real Ue,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real (&phiqq)[Nsound]
);

void calc_WPS_TNO(
    const Real delta,
    Real tauWall,
    const Real edgeVel,
    const Real (&omega)[Nsound],
    const Real rho,
    const Real nu,
    const int isSuction,
    Real (&phiqq)[Nsound]);

//Real calc_Spp_Freq(
//    Real c0,  Real rho0, Real C,Real MX, Real omega,
//    Real X, Real Y, Real Z,
//    Real S, Real Phi_qq_input,int Order
//);


void TE_noise_outer(
    // Flow / geometry parameters formerly in 'inputs':
    Real M, Real U, Real x, Real y, Real z,
    Real b, Real Ky, Real c, Real span,

    // Fluid properties formerly in 'fluid':
    Real c0, Real rho, Real nu,

    // Array of frequencies:
    const Real omega[Nsound],
    const Real Ue_bot, const Real Ue_top,
    Real (&WPS_lower)[Nsound], Real (&WPS_upper)[Nsound],

    // Output wall-pressure spectrum:
    Real (&farfieldSpectra)[Nsound]
);

#endif
