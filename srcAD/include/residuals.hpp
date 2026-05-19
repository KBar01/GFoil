#pragma once

// get_funcs.hpp and vector_ops.hpp must be included before residuals_shared.hpp.
#include "data_structs.hpp"
#include "get_funcs.hpp"
#include "vector_ops.hpp"
// Shared physics (upwind, upwind_half, residual_station … residual_transition_forced):
#include "residuals_shared.hpp"

// Isolc and Isolv are genuine template structs; the others are now template
// aliases defined in data_structs.hpp via data_structs_shared.hpp.
template<typename Real> struct Isolc;
template<typename Real> struct Isolv;

// AD-specific: wake_sys, residual_transition, residual_transition_forced
// use Vsol<Real>/Foil<Real>/Glob<Real> (template structs from data_structs.hpp)
// which differ from the non-template equivalents in src/include/data_structs.h.

