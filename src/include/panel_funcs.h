// Thin forwarding header for the fwd build.
// Template panel functions (panel_info, panel_linvortex_stream, etc.) are in panel_funcs.hpp.
// inviscid_velocity and dvelocity_dgamma are now template<Real,FoilT> in panel_funcs.hpp;
// the declarations below are the non-template fwd wrappers implemented in panel_funcs.cpp.
#pragma once

#include "real_type.h"
#include "panel_funcs.hpp"

struct Foil;  // full definition in data_structs.h; only a reference needed here

void inviscid_velocity(const Foil& foil, const Real* gammas, const Real& Vinf,
    const Real& alpha, const Real& CPx, const Real& CPy, Real* velocity);

void dvelocity_dgamma(const Foil& foil, const Real& CPx, const Real& CPy, Real* V_G);
