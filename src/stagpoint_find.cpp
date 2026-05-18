#include <iostream>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "solver_funcs.hpp"


void stagpoint_find(Isol& isol, const Foil& foil, const Wake& wake) {

    int j=0;
    // Find first positive gamma
    for (; j < Ncoords; ++j) {
        if (isol.gammas[j] > 0) break;
    }
    int I[2] = {j-1, j};
    isol.stagIndex[0] = I[0];
    isol.stagIndex[1] = I[1];

    Real G[2] = {isol.gammas[I[0]], isol.gammas[I[1]]};
    Real S[2] = {foil.s[I[0]], foil.s[I[1]]};

    Real den = G[1] - G[0];
    Real w1 = G[1] / den;
    Real w2 = -G[0] / den;

    isol.stagArcLocation = w1 * S[0] + w2 * S[1];

    isol.stagXLocation[0] = foil.x[colMajorIndex(0,j-1,2)] * w1 + foil.x[colMajorIndex(1,j,2)] * w2;
    isol.stagXLocation[1] = foil.x[colMajorIndex(1,j-1,2)] * w1 + foil.x[colMajorIndex(1,j,2)] * w2;

    Real st_g1 = G[1] * (S[0] - S[1]) / (den * den);
    isol.sstag_g[0] = st_g1;
    isol.sstag_g[1] = -st_g1;

    for (int i = j; i < Ncoords; ++i) {
        isol.edgeVelSign[i] = 1.0;
    }

    for (int i = 0; i < Ncoords; ++i) {
        isol.distFromStag[i] = std::abs(foil.s[i] - isol.stagArcLocation);
    }

    for (int i = 0; i < Nwake; ++i) {
        isol.distFromStag[Ncoords + i] = wake.s[i] - isol.stagArcLocation;
    }
}

// range() is now inline in solver_funcs.hpp; this is the fwd non-template wrapper.
void identify_surfaces(const Isol& isol, Vsol& vsol) {
    identify_surfaces<>(isol, vsol);
}
