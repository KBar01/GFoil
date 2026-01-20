#include <iostream>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"


void stagpoint_find(Isol& isol, const Foil& foil,const Wake& wake) {


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

    // Compute x location in column-major form
    isol.stagXLocation[0] = foil.x[colMajorIndex(0,j-1,2)] * w1 + foil.x[colMajorIndex(1,j,2)] * w2;   // x-coordinate
    isol.stagXLocation[1] = foil.x[colMajorIndex(1,j-1,2)] * w1 + foil.x[colMajorIndex(1,j,2)] * w2; // y-coordinate

    // Compute linearization
    Real st_g1 = G[1] * (S[0] - S[1]) / (den * den);
    isol.sstag_g[0] = st_g1;
    isol.sstag_g[1] = -st_g1;

    // Compute sign conversion (sgnue)
    for (int i = j; i < Ncoords; ++i) {
        isol.edgeVelSign[i] = 1.0;
    }

    // Compute distance from stagnation point (xi)
    for (int i = 0; i < Ncoords; ++i) {
        isol.distFromStag[i] = std::abs(foil.s[i] - isol.stagArcLocation);
    }

    for (int i = 0; i < Nwake; ++i) {
        isol.distFromStag[Ncoords + i] = wake.s[i] - isol.stagArcLocation;
    }

}

// Helper function to mimic Python's range(start, end, step)
std::vector<int> range(int start, int end, int step = 1) {
    std::vector<int> result;
    if (step > 0) {
        for (int i = start; i < end; i += step)
            result.push_back(i);
    } else if (step < 0) {
        for (int i = start; i > end; i += step)
            result.push_back(i);
    }
    return result;
}


// Core function to fill in surface indices
void identify_surfaces(const Isol& isol, Vsol& vsol) {
    vsol.Is.clear(); // Clear any previous data

    // Lower surface (reverse order from Istag[0] to 0)
    vsol.Is.push_back(range(isol.stagIndex[0], -1, -1));

    // Upper surface (from Istag[1] to foil.N-1)
    vsol.Is.push_back(range(isol.stagIndex[1], Ncoords));

    // Wake (from foil.N to foil.N + wake.N - 1)
    vsol.Is.push_back(range(Ncoords, Ncoords+Nwake));
}
