// Forwarding header for the forward solver.
// residual_station and residual_station_forced are templatised in residuals.hpp;
// the non-template wrappers below allow existing call sites to pass int literals
// for Real parameters without explicit template instantiation.
#pragma once

#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "get_funcs.h"     // must precede residuals_shared.hpp so get_H etc. are visible at instantiation
#include "vector_ops.hpp"  // must precede residuals_shared.hpp so cnp:: ops are available
#include "residuals_shared.hpp"

// Non-template declarations — struct-dependent functions kept in src/residuals.cpp:

void wake_sys(const Vsol& vsol, const Foil& foil, const Glob& glob, const Param& param,
    Real (&R)[3], Real (&R_U)[36], int (&J)[3]);

void residual_transition(
    const Real* U1, const Real* U2,
    const Real x1, const Real x2,
    const Real aux1, const Real aux2,
    const Param& param,
    Real (&R)[3], Real (&R_U)[24], Real (&R_x)[6]);

void residual_transition_forced(
    const Real* U1, const Real* U2,
    const Real x1, const Real x2,
    const Param& param,
    const Real transPos,
    Real (&R)[3], Real (&R_U)[24], Real (&R_x)[6]);

// Non-template wrappers: allow int→Real implicit conversion at existing call sites
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

inline void residual_station_forced(
    const Real* U1, const Real* U2,
    const Real x1, const Real x2,
    const Real aux1, const Real aux2,
    const bool wake, const bool turb, const bool simi,
    const Param& param,
    const Real& ncrit,
    Real (&R)[3], Real (&R_U)[24], Real (&R_x)[6])
{
    ::residual_station_forced<true>(U1, U2, x1, x2, aux1, aux2,
                                    wake, turb, simi, param, ncrit, R, R_U, R_x);
}
