#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.hpp"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.hpp"
#include "vector_ops.hpp"
#include "main_func.h"
#include "solver_funcs.hpp"
#include "restart_state.h"
#include <chrono>
#include <fstream>

#include "nlohmann/json.hpp"

using json = nlohmann::json;


Real euc_norm(const Real* R, int size) {
    Real sum = 0.0;
    for (int i = 0; i < size; ++i) {
        sum += R[i] * R[i];
    }
    return std::sqrt(sum);
}

bool solve_coupled(const Oper& oper, const Foil& foil, const Wake& wake,
    Param& param, Vsol& vsol, Isol& isol, Glob& glob,
    RestartState* restartOut,
    std::string* failure_mode_out) {

    int nNewton = param.niglob;
    bool converged = false;
    constexpr int Rsize = 3*(Ncoords + Nwake);
    constexpr int Rallsize = 4*(Ncoords + Nwake);

    // GFOIL_DEBUG=1 enables per-iteration residual/omega/transition diagnostics.
    bool debugMode = (std::getenv("GFOIL_DEBUG") != nullptr);

    // Return last-laminar node index (0-based in Is[si]) for surface si.
    // Returns Is.size()-1 when fully laminar, -1 when fully turbulent from node 0.
    auto find_ilam = [&](int si) -> int {
        if (si >= static_cast<int>(vsol.Is.size())) return -1;
        const auto& Is = vsol.Is[si];
        for (int k = 0; k < static_cast<int>(Is.size()); ++k) {
            if (vsol.turb[Is[k]]) return k - 1;
        }
        return static_cast<int>(Is.size()) - 1;
    };

    // Return amp/ctau value at surface si, node index k (index into Is[si]).
    auto get_amp = [&](int si, int k) -> double {
        if (si >= static_cast<int>(vsol.Is.size())) return 0.0;
        const auto& Is = vsol.Is[si];
        if (k < 0 || k >= static_cast<int>(Is.size())) return 0.0;
        return glob.U[colMajorIndex(2, Is[k], 4)].getValue();
    };

    if (debugMode) {
        std::fprintf(stderr,
            "DBG  iter    L2_resid      omega  ilam_bot ilam_top  amp_bot  amp_top\n");
    }

    // Failure-mode tracking: detect period-2 oscillation at the transition front.
    int ilam_bot_prev = -1, ilam_top_prev = -1;
    int stable_iters = 0;
    int oscillation_count = 0;
    double resid_prev = 1e20;

    for (int i = 0; i < 60; ++i) {

        // Main loop solving coupled system

        build_glob_RV(foil, vsol, isol, glob, param);
        Real residualNorm = euc_norm(glob.R, Rsize);


        if (residualNorm < param.rtol) {

            if (debugMode) {
                int ib = find_ilam(0), it = find_ilam(1);
                std::fprintf(stderr,
                    "DBG  %4d  %12.5e  CONVERGED  ilam_bot=%d ilam_top=%d\n",
                    i, residualNorm.getValue(), ib, it);
            }

            solve_glob(foil,isol,glob,vsol,oper,0);

            std::vector<double> states_vec(RVdimension);
            for (int k = 0; k < RVdimension; ++k)
                states_vec[k] = glob.U[k].getValue();

            std::vector<double> jac_vec(glob.R_V_latest);
            std::vector<int> jac_row_vec(glob.R_V_latest);
            std::vector<int> jac_col_vec(glob.R_V_latest);
            for (int k = 0; k < glob.R_V_latest; ++k) {
                jac_vec[k]     = glob.R_V_vals[k].getValue();
                jac_row_vec[k] = glob.R_V_rows[k];
                jac_col_vec[k] = glob.R_V_cols[k];
            }

            if (restartOut != nullptr) {
                restartOut->states = states_vec;
                restartOut->turb.assign(vsol.turb, vsol.turb + Ncoords + Nwake);
                restartOut->stag   = {isol.stagIndex[0], isol.stagIndex[1]};
                restartOut->RVnz   = glob.R_V_latest;
                restartOut->RVvals = jac_vec;
                restartOut->RVrows = jac_row_vec;
                restartOut->RVcols = jac_col_vec;
            } else {
                json restart;
                restart["states"] = states_vec;
                restart["turb"]   = vsol.turb;
                restart["stag"]   = isol.stagIndex;
                restart["RVvals"] = jac_vec;
                restart["RVrows"] = jac_row_vec;
                restart["RVcols"] = jac_col_vec;
                restart["RVnz"]   = glob.R_V_latest;
                std::ofstream fout("restart.json");
                fout << restart.dump(4);
            }

            clear_RV(glob, isol, vsol, foil, param);
            converged = true;
            glob.convergenceIteration = i;
            break;
        }
        
        solve_glob(foil, isol, glob, vsol, oper, 1);
        Real omega = update_state(oper, param, glob, vsol);
        clear_RV(glob, isol, vsol, foil, param);
        for (int entry = 0; entry < Rallsize; ++entry) {
            glob.R[entry] = 0;
        }
        stagpoint_move(isol, glob, foil, wake, vsol);
        update_transition(glob, vsol, isol, param, i);

        {
            int ib = find_ilam(0), it = find_ilam(1);
            if (debugMode) {
                double ab = (ib >= 0) ? get_amp(0, ib) : 0.0;
                double at = (it >= 0) ? get_amp(1, it) : 0.0;
                std::fprintf(stderr,
                    "DBG  %4d  %12.5e  %8.5f  %8d %8d  %8.4f  %8.4f\n",
                    i, residualNorm.getValue(), omega.getValue(), ib, it, ab, at);
            }

            bool ilam_stable = (ib == ilam_bot_prev && it == ilam_top_prev);
            if (ilam_stable) {
                ++stable_iters;
                if (residualNorm.getValue() > resid_prev) ++oscillation_count;
            } else {
                stable_iters    = 0;
                oscillation_count = 0;
            }
            ilam_bot_prev = ib;
            ilam_top_prev = it;
            resid_prev    = residualNorm.getValue();
        }
    }

    if (!converged && failure_mode_out != nullptr) {
        double rn = resid_prev;  // residual at the last iteration
        if (stable_iters >= 15 && oscillation_count >= 6 && rn < 1.0)
            *failure_mode_out = "transition_front_oscillation";
        else if (rn >= 1.0)
            *failure_mode_out = "diverged";
        else
            *failure_mode_out = "no_convergence";
    }

    return converged;
}