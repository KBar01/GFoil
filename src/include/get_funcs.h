// Forwarding header for the forward solver.
// All BL physics functions (get_uk, get_H, get_cf, etc.) are now
// template<typename Real, typename ParamT> in get_funcs.hpp.
// get_ueinv is kept as a non-template because it takes Isol (monolithic),
// while srcAD uses Isolc<Real>+Isolv<Real> — the Isol/Isolc split prevents merging.
#pragma once

#include "real_type.h"
#include "data_structs.h"
#include "get_funcs.hpp"

// Kept separate from get_funcs.hpp: requires Isol (forward) vs Isolc+Isolv (AD)
void get_ueinv(const Isol& isol, Real* ueinv);
