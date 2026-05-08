// All BL physics functions are now in src/include/get_funcs.hpp (shared templates).
// get_ueinv is kept here because it takes the non-template Isol struct (forward solver),
// while the AD solver uses Isolc<Real> + Isolv<Real> for correct tape segmentation.
// Merging would require resolving the Isol/Isolc split in data_structs.
#include "real_type.h"
#include "data_structs.h"
#include "get_funcs.h"

void get_ueinv(const Isol& isol, Real* ueinv) {
    for (int i = 0; i < Ncoords; ++i) {
        ueinv[i] = isol.edgeVelSign[i] * isol.gammas[i];
    }
    for (int i = 0; i < Nwake; ++i) {
        ueinv[Ncoords+i] = isol.uewi[i];
    }
    ueinv[Ncoords] = ueinv[Ncoords-1]; // continuity of upper surface and wake velocity
}
