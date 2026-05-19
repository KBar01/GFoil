#pragma once

// get_funcs.hpp and vector_ops.hpp must be included before residuals_shared.hpp.
#include "data_structs.hpp"
#include "get_funcs.hpp"
#include "vector_ops.hpp"
// Shared physics (upwind, upwind_half, residual_station, residual_transition, wake_sys):
#include "residuals_shared.hpp"

// Isolc and Isolv are genuine template structs; the others are now template
// aliases defined in data_structs.hpp via data_structs_shared.hpp.
template<typename Real> struct Isolc;
template<typename Real> struct Isolv;

