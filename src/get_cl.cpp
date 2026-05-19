#include <iostream>
#include <cmath>
#include "real_type.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "solver_funcs.hpp"

void calc_force(const Oper& op, const Geom& geom, const Param& par,
                const Isol& isol, const Foil& foil,
                const Glob& glob, Post& post) {
    // isol is unused in this function — it was always dead weight in the signature
    calc_force<>(op, geom, par, foil, glob, post);
}
