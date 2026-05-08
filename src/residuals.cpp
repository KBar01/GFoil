#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"
#include "residuals.h"  // pulls in residuals.hpp + non-template wrappers for internal calls


void wake_sys(const Vsol& vsol, const Foil& foil, const Glob& glob, const Param& param,
    Real (&R)[3], Real (&R_U)[36], int (&J)[3])
{
    wake_sys<true, Real>(vsol, foil, glob, param, R, R_U, J);
}

void residual_transition(
    const Real* U1, const Real* U2,
    const Real x1, const Real x2,
    const Real aux1, const Real aux2,
    const Param& param,
    Real (&R)[3], Real (&R_U)[24], Real (&R_x)[6])
{
    residual_transition<true, Real>(U1, U2, x1, x2, aux1, aux2, param, R, R_U, R_x);
}


void residual_transition_forced(
    const Real* U1, const Real* U2,
    const Real x1, const Real x2,
    const Param& param,
    const Real transPos,
    Real (&R)[3], Real (&R_U)[24], Real (&R_x)[6])
{
    residual_transition_forced<true, Real>(U1, U2, x1, x2, param, transPos, R, R_U, R_x);
}