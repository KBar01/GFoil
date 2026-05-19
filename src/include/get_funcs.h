// Forwarding header for the forward solver.
// All BL physics functions (get_uk, get_H, get_cf, etc.) are now
// template<typename Real, typename ParamT> in get_funcs.hpp.
// get_ueinv is kept as a non-template because it takes Isol (monolithic),
// while srcAD uses Isolc<Real>+Isolv<Real> — the Isol/Isolc split prevents merging.
#pragma once

#include "real_type.h"
#include "data_structs.h"
#include "get_funcs.hpp"

// get_ueinv: fwd non-template inline (Isol is a unified struct; AD uses isolc+isolv).
// Guard ensures Isol/Real are defined (they come from data_structs.h above).
#ifdef AIRFOIL_STRUCTS_H
inline void get_ueinv(const Isol& isol, Real* ueinv) {
    for (int i = 0; i < Ncoords; ++i)
        ueinv[i] = isol.edgeVelSign[i] * isol.gammas[i];
    for (int i = 0; i < Nwake; ++i)
        ueinv[Ncoords+i] = isol.uewi[i];
    ueinv[Ncoords] = ueinv[Ncoords-1];
}
#endif // AIRFOIL_STRUCTS_H
