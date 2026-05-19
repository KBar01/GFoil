#include <codi.hpp>
#include <Eigen/Dense>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.hpp"
#include "data_structs.h"
#include "solve_inv.hpp"

// colMajorIndex: canonical non-inline definition for the forward solver.
// Satisfies the declaration in main_func.h used by all other src/ TUs.
// (An inline version also exists in sparselinsolve.hpp for files that include that header.)
int colMajorIndex(int row, int col, int num_rows) {
    return row + col * num_rows;
}

// Non-template wrapper: allows existing call sites in the forward solver
// (e.g. src/main.cpp) to call build_gamma_codi without explicit template
// arguments, and satisfies the declaration in src/include/main_func.h.
// The actual physics lives in the shared template in solve_inv.hpp.
void build_gamma_codi(Isol& isol, const Foil& foil, const Oper& op)
{
    build_gamma_codi<Real>(isol, foil, op);
}
