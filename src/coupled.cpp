#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"
#include "main_func.h"
#include <chrono>
#include <fstream>

#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;


Real euc_norm(const Real* R, int size) {
    Real sum = 0.0;
    for (int i = 0; i < size; ++i) {
        sum += R[i] * R[i];
    }
    return std::sqrt(sum);
}

bool solve_coupled(const Oper& oper, const Foil& foil, const Wake& wake,
    Param& param, Vsol& vsol, Isol& isol, Glob& glob, Trans&tdata, const bool force) {

    int nNewton = param.niglob;
    bool converged = false;
    constexpr int Rsize = 3*(Ncoords + Nwake);
    constexpr int Rallsize = 4*(Ncoords + Nwake);

    for (int i = 0; i < 60; ++i) {
        
        // Main loop solving coupled system


        build_glob_RV(foil, vsol, isol, glob, param,tdata);
        Real residualNorm = euc_norm(glob.R, Rsize);


        if (residualNorm < param.rtol) {
            
            // if converged return states and also dR/dU info for AD code
            solve_glob(foil,isol,glob,vsol,oper,0);
            json restart;
            std::vector<double> states_vec(RVdimension);
            for (int k = 0; k < RVdimension; ++k){
                states_vec[k] = glob.U[k].getValue();
            }

            restart["states"] = states_vec;
            restart["turb"]   = vsol.turb;
            restart["stag"] = isol.stagIndex;

            std::vector<double> jac_vec(glob.R_V_latest);
            std::vector<int> jac_row_vec(glob.R_V_latest);
            std::vector<int> jac_col_vec(glob.R_V_latest);
            for (int k = 0; k < glob.R_V_latest; ++k){
                jac_vec[k] = glob.R_V_vals[k].getValue();
                jac_row_vec[k] = glob.R_V_rows[k] ;
                jac_col_vec[k] = glob.R_V_cols[k] ;
            }
            restart["RVvals"] = jac_vec;
            restart["RVrows"] = jac_row_vec;
            restart["RVcols"] = jac_col_vec;
            restart["RVnz"] = glob.R_V_latest;
            // Write JSON file
            std::ofstream fout("restart.json");
            fout << restart.dump(4); 
            
            clear_RV(glob, isol, vsol, foil, param);
            converged = true;
            glob.convergenceIteration = i;
            break;
        }
        
        solve_glob(foil, isol, glob, vsol, oper, 1);
        update_state(oper, param, glob, vsol);
        clear_RV(glob, isol, vsol, foil, param);
        for (int entry = 0; entry < Rallsize; ++entry) {
            glob.R[entry] = 0;
        }
        stagpoint_move(isol, glob, foil, wake, vsol);
        update_transition(glob, vsol, isol, param, tdata, force);

    }

    return converged;
}

//#endif