#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "vector_ops.hpp"
#include "main_func.h"
#include "data_structs.h"
#include "solver_funcs.hpp"


void stagpoint_move(Isol& isol, Glob& glob, const Foil& foil, const Wake& wake, Vsol& vsol) {

    // Preamble: determine which stagnation panel I[0..1] we are on.
    int* I = isol.stagIndex;
    if (glob.U[colMajorIndex(3,I[1],4)] < 0) {

        int j;
        for (j = I[1]; j < Ncoords; ++j) {
            if (glob.U[colMajorIndex(3,j,4)] > 0) break;
        }
        int I1 = j;
        for (j = I[1]; j < I1; ++j) { glob.U[colMajorIndex(3,j,4)] *= -1.0; }
        I[0] = I1 - 1;
        I[1] = I1;
    }
    else if (glob.U[colMajorIndex(3,I[0],4)] < 0) {

        int j;
        for (j = I[0]; j >= 0; --j) {
            if (glob.U[colMajorIndex(3,j,4)] > 0) break;
        }
        assert(j > 0 && "no stagnation point");
        int I0 = j;
        for (j = I0+1; j <= I[0]; ++j) { glob.U[colMajorIndex(3,j,4)] *= -1.0; }
        I[0] = I0;
        I[1] = I0 + 1;
    }

    // Core: update all stagnation-dependent state (always safe to call).
    int stagPanel[2] = {I[0], I[1]};
    stagpoint_move_impl<Real>(isol, glob, foil, wake, vsol, stagPanel);
}
