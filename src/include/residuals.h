#pragma once

#include "real_type.h"
#include "panel_funcs.hpp"
#include "data_structs.h"
#include "get_funcs.hpp"     // must precede residuals_shared.hpp so get_H etc. are visible at instantiation
#include "vector_ops.hpp"  // must precede residuals_shared.hpp so cnp:: ops are available
#include "residuals_shared.hpp"

// Non-template wrapper: allow int→Real implicit conversion at existing call sites
// (template deduction fails when int literals are passed for Real parameters).

inline void residual_station(
    const Real* U1, const Real* U2,
    const Real x1, const Real x2,
    const Real aux1, const Real aux2,
    const bool wake, const bool turb, const bool simi,
    const Param& param,
    Real (&R)[3], Real (&R_U)[24], Real (&R_x)[6])
{
    ::residual_station<true>(U1, U2, x1, x2, aux1, aux2,
                             wake, turb, simi, param, R, R_U, R_x);
}
// residual_station_forced removed: forced transition deleted.
