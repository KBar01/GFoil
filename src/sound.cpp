// src/sound.cpp — forward-solver acoustic post-processing.
//
// The shared template implementations of calc_WPS and calc_OASPL live in
// src/include/sound.hpp.  This file provides a non-template wrapper that
// preserves the existing call signature in src/main.cpp (WPSjson before model)
// and routes to the WriteJSON=true instantiation.

#include <cmath>
#include "data_structs.h"
#include "real_type.h"
#include "sound.hpp"

// Non-template wrapper: backward-compatible signature for src/main.cpp.
// WPSjson is passed before model here (matching the original API) and
// forwarded to the WriteJSON=true template instantiation.
Real calc_OASPL(const Real* botStates, const Real* topStates,
                const Real chordScale, const Real Uinf,
                const Real X, const Real Y, const Real Z,
                const Real S, const Real nu, const Real rho,
                const int WPSjson, const std::string& model)
{
    return calc_OASPL<Real, true>(botStates, topStates, chordScale, Uinf,
                                   X, Y, Z, S, nu, rho, model, WPSjson);
}
