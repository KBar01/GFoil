#include <iostream>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "solver_funcs.hpp"

void space_wake_nodes(const Real& wakeLength, const Real& firstPanelLength,
                      Real* wakeSpacing, const Foil& foil, Wake& wake) {
    space_wake_nodes<>(wakeLength, firstPanelLength, wakeSpacing, foil, wake);
}



void build_wake(const Foil& foil, const Geom& geom, const Oper& op, Isol& isol, Wake& wake) {
    build_wake_impl<>(foil, geom, op, isol, wake);
}


void set_wake_gap(const Foil& foil, const Isol& isol, Vsol& vsol) {
    set_wake_gap<>(foil, isol, vsol);
}