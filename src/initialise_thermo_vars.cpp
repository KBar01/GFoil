#include <iostream>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "solver_funcs.hpp"

void init_thermo(const Oper& oper, Param& param, const Geom& geom) {
    init_thermo<>(oper, param, geom);
}
