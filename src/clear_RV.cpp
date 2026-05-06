#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"


void clear_RV(Glob&glob, const Isol&isol,const Vsol&vsol, const Foil&foil,const Param&param){

    glob.R_V_vals.clear();
    glob.R_V_rows.clear();
    glob.R_V_cols.clear();

}

