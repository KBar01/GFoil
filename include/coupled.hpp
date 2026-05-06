#pragma once
#include <fstream>
#include "real_type.hpp"
#include "data_structs.hpp"
#include "build_global_sys.hpp"
#include "solve_glob.hpp"
#include "clear_RV.hpp"
#include "stagmove.hpp"
#include "update_state.hpp"
#include "update_transition.hpp"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

template<typename Real>
Real euc_norm(const Real* R, int size) {
    Real sum = 0.0;
    for (int i=0;i<size;++i) sum += R[i]*R[i];
    return std::sqrt(sum);
}

template<typename Real>
bool solve_coupled(const Oper<Real>& oper, const Foil<Real>& foil, const Wake<Real>& wake,
    Param<Real>& param, Vsol<Real>& vsol, const Isolc<Real>& isolc, Isolv<Real>& isol,
    Glob<Real>& glob, Trans<Real>& tdata, const bool force)
{
    bool converged = false;
    const int Rsize    = 3*(glob.nc + glob.nw);
    const int Rallsize = 4*(glob.nc + glob.nw);

    for (int i=0; i<60; ++i) {
        build_glob_RV<Real>(foil, vsol, isol, glob, param, tdata);
        Real residualNorm = euc_norm<Real>(glob.R.data(), Rsize);

        if (residualNorm < param.rtol) {
            solve_glob<Real>(foil, isolc, isol, glob, vsol, oper, 0);

            json restart;
            std::vector<double> states_vec(RVdimension);
            for (int k=0;k<RVdimension;++k) states_vec[k] = glob.U[k].getValue();
            restart["states"] = states_vec;
            restart["turb"]   = vsol.turb;
            restart["stag"]   = isol.stagIndex;

            const int nnz = (int)glob.R_V_vals.size();
            std::vector<double> jac_vec(nnz);
            std::vector<int>    jac_row(nnz), jac_col(nnz);
            for (int k=0;k<nnz;++k) {
                jac_vec[k] = glob.R_V_vals[k].getValue();
                jac_row[k] = glob.R_V_rows[k];
                jac_col[k] = glob.R_V_cols[k];
            }
            restart["RVvals"] = jac_vec;
            restart["RVrows"] = jac_row;
            restart["RVcols"] = jac_col;
            restart["RVnz"]   = nnz;

            std::ofstream fout("restart.json");
            fout << restart.dump(4);

            clear_RV<Real>(glob);
            converged = true;
            glob.convergenceIteration = i;
            break;
        }

        solve_glob<Real>(foil, isolc, isol, glob, vsol, oper, 1);
        update_state<Real>(oper, param, glob, vsol);
        clear_RV<Real>(glob);
        for (int entry=0;entry<Rallsize;++entry) glob.R[entry] = 0;
        stagpoint_move<Real>(isol, glob, foil, wake, vsol);
        update_transition<Real>(glob, vsol, isol, param, tdata, force);
    }

    return converged;
}
