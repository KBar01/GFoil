#pragma once
#include "real_type.h"
#include "data_structs.h"
#include "restart_state.h"
#include <string>

// When fwdOut != nullptr: skip out.json write, fill struct instead.
// When restartOut != nullptr: skip restart.json write, fill struct instead.
// When warmStart != nullptr: initialise glob.U/vsol.turb from in-memory state
//   instead of reading restart.json (pybind11 warm-start path).
bool runCode(
    bool fromRestart,
    const Real nCrit,
    const Real Ufac,
    const Real TEfac,
    const Real chordScaling,
    const Real (&inXcoords)[Nin],
    Real (&inYcoords)[Nin],
    Real alphad,
    Real Re,
    Real Ma,
    Real rhoInf,
    Real kinViscInf,
    const std::string model,
    const Real sampleTE,
    const Real* obsX,
    const Real* obsY,
    const Real* obsZ,
    int nObs,
    const Real S,
    RestartState* restartOut             = nullptr,
    ForwardResult* fwdOut                = nullptr,
    const RestartState* warmStart        = nullptr,
    int aWeighting                       = 0,
    Real ncrithyst                       = 0.2,
    bool verbose                         = false,
    double f_min                         = 200.0,
    double f_max                         = 20000.0,
    double xft_lower                     = 1.0,
    double xft_upper                     = 1.0
);
