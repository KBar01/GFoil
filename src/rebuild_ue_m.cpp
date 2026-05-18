#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "solver_funcs.hpp"

void rebuild_ue_m(const Foil& foil, const Wake& wake,
                  const Isol& isol, Vsol& vsol, bool realloc) {
    rebuild_ue_m<>(foil, wake, isol, vsol, realloc);
}
