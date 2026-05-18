#include <iostream>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "solver_funcs.hpp"


void stagpoint_find(Isol& isol, const Foil& foil, const Wake& wake) {
    // isol serves as both gamma source and variable destination (unified fwd struct)
    stagpoint_find_impl<true>(isol, isol, foil, wake);
}

// range() is now inline in solver_funcs.hpp; this is the fwd non-template wrapper.
void identify_surfaces(const Isol& isol, Vsol& vsol) {
    identify_surfaces<>(isol, vsol);
}
