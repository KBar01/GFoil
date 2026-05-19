// Non-template wrapper functions for the fwd build.
// Template implementations (FoilT duck-typed) now live in panel_funcs.hpp.

#include "real_type.h"
#include "data_structs.h"
#include "panel_funcs.hpp"

void inviscid_velocity(const Foil& foil, const Real* gammas, const Real& Vinf,
    const Real& alpha, const Real& CPx, const Real& CPy, Real* velocity)
{
    inviscid_velocity<>(foil, gammas, Vinf, alpha, CPx, CPy, velocity);
}

void dvelocity_dgamma(const Foil& foil, const Real& CPx, const Real& CPy, Real* V_G)
{
    dvelocity_dgamma<>(foil, CPx, CPy, V_G);
}
