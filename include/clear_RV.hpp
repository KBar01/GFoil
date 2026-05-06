#pragma once
#include "data_structs.hpp"

template<typename Real>
void clear_RV(Glob<Real>& glob) {
    glob.R_V_vals.clear();
    glob.R_V_rows.clear();
    glob.R_V_cols.clear();
}
