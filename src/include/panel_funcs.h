// Thin forwarding header — template panel functions are defined in panel_funcs.hpp.
// This header additionally declares the non-template functions that depend on
// the Foil struct (inviscid_velocity, dvelocity_dgamma), which are implemented
// in src/panel_funcs.cpp and cannot be shared as templates.
#pragma once

#include "real_type.h"
#include "panel_funcs.hpp"

struct Foil;  // full definition in data_structs.h; only a reference needed here

void inviscid_velocity(const Foil& foil, const Real* gammas, const Real& Vinf,
    const Real& alpha, const Real& CPx, const Real& CPy, Real* velocity);

void dvelocity_dgamma(const Foil& foil, const Real& CPx, const Real& CPy, Real* V_G);
