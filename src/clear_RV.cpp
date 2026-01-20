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

    int lastIndex = glob.R_V_latest;

    for (int i=0;i<lastIndex;++i){
        glob.R_V_vals[i] = 0.0;
        glob.R_V_rows[i] = 0;
        glob.R_V_cols[i] = 0;
    }

    glob.R_V_latest = 0;

}

