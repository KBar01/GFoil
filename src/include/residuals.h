#ifndef RESIDUALS_H
#define RESIDUALS_H


#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"



void wake_sys(const Vsol& vsol, const Foil& foil, const Glob& glob, const Param& param,
    Real (&R)[3], Real (&R_U)[36], int (&J)[3]);


Real upwind_half(

    // Sets upw to 0.5 and upw_U to 0 
    const Real f1,
    const Real (&f1_U1)[4],
    const Real f2,
    const Real (&f2_U2)[4],
    Real (&f_U)[8]);

void residual_station(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Real aux1,
    const Real aux2,
    const bool wake,
    const bool turb,
    const bool simi,
    const Param&param,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6]);


void residual_transition(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Real aux1,
    const Real aux2,
    const Param& param,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6]
);


void residual_transition_forced(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Param& param,
    const Real transPos,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6]
);


#endif
