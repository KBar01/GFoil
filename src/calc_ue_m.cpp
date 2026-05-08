// calc_ue_m, compute_Dw, and solve_sys_ue are now template<typename Real, ...>
// in src/include/calc_ue_m.hpp.  This file provides the non-template wrapper
// that satisfies the declaration in src/include/main_func.h and allows existing
// call sites to omit explicit template arguments.
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "vector_ops.hpp"
#include "calc_ue_m.hpp"

void calc_ue_m(const Foil& foil, const Wake& wake, Isol& isol, Vsol& vsol)
{
    calc_ue_m<Real>(foil, wake, isol, vsol);
}
