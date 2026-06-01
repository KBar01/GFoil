#include <iostream>
#include <cmath>
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

    bool converged = false;
    constexpr int Rsize = 3*(Ncoords + Nwake);
    constexpr int Rallsize = 4*(Ncoords + Nwake);

    // Return last-laminar node index (0-based in Is[si]) for surface si.
    auto find_ilam = [&](int si) -> int {
        if (si >= static_cast<int>(vsol.Is.size())) return -1;
        const auto& Is = vsol.Is[si];
        for (int k = 0; k < static_cast<int>(Is.size()); ++k) {
            if (vsol.turb[Is[k]]) return k - 1;
        }
        return static_cast<int>(Is.size()) - 1;
    };

    // Per-surface transition tracking and ctau-freeze cycle detection.
    // All plain doubles/ints/bools — no Real — so no CoDi tape contamination.
    double resid_buf[8]         = {};        // circular buffer of last 8 L2 norms
    int    resid_pos            = 0;         // total entries written (head = pos % 8)
    int    ilam_prev[2]         = {-1, -1};  // ilam from previous iter per surface
    int    stable_ilam_iters[2] = {0, 0};    // consecutive iters with ilam unchanged
    bool   ctau_freeze[2]       = {false, false}; // freeze flag per surface
    double prev_amp[2]          = {-1.0, -1.0};   // last amp at Is[ilam0] per surface
    // prev_ctau[si][k]: last Newton-updated ctau at Is[ilam0+1+k] for surface si.
    // Averaging the first 3 turbulent nodes damps the multi-node oscillation.
    double prev_ctau[2][3]      = {{-1.0,-1.0,-1.0},{-1.0,-1.0,-1.0}};

    int    nan_skip_count = 0;    // consecutive iters matching the NaN-lock signature
    double prev_resid_val = -1.0; // residual value from the previous iteration
    bool   early_exit     = false; // set true when nan_lock forces early termination

    for (int i = 0; i < 60; ++i) {

        build_glob_RV(foil, vsol, isol, glob, param);
        Real residualNorm = euc_norm(glob.R, Rsize);


        if (residualNorm < param.rtol) {

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
        update_transition(glob, vsol, isol, param, foil, i, ctau_freeze[1], ctau_freeze[0],
                          prev_amp + 1, prev_amp, prev_ctau[1], prev_ctau[0]);

        {
            int cur_ilam[2] = {find_ilam(0), find_ilam(1)};

            // Update residual circular buffer (plain double, no CoDi).
            int slot = resid_pos % 8;
            double oldest_resid = resid_buf[slot];
            resid_buf[slot] = residualNorm.getValue();
            bool buf_full = (resid_pos >= 8);
            ++resid_pos;

            // Cycling: current residual has not improved by 2x vs 8 iterations ago.
            bool not_improving = buf_full &&
                (residualNorm.getValue() > oldest_resid / 2.0);

            // Per-surface: update stable counter and freeze flags.
            for (int si = 0; si < 2; ++si) {
                if (cur_ilam[si] != ilam_prev[si]) {
                    stable_ilam_iters[si] = 0;
                    ctau_freeze[si] = false;      // reset when ilam moves
                    prev_amp[si]  = -1.0;
                    prev_ctau[si][0] = prev_ctau[si][1] = prev_ctau[si][2] = -1.0;
                } else {
                    ++stable_ilam_iters[si];
                }
                ilam_prev[si] = cur_ilam[si];

                // Activate freeze when transition has been stable ≥8 iters,
                // residual is cycling (not improving), and not yet converged.
                if (stable_ilam_iters[si] >= 8 && not_improving &&
                    residualNorm.getValue() > param.rtol.getValue()) {
                    ctau_freeze[si] = true;
                }
                // Keep freeze active once triggered — releasing early resets
                // the averaged history and restarts the same oscillation cycle.
                // The only release paths are: ilam moves (transition shifts) or
                // the solver converges (residualNorm < rtol → loop exits).
            }

            // Co-activate: once either surface triggers a freeze, also freeze the
            // partner surface if it has had ≥4 consecutive stable-ilam iterations.
            // This prevents one surface cycling freely while the other is frozen,
            // which creates cross-coupling oscillations in the global Newton system.
            if (ctau_freeze[0] || ctau_freeze[1]) {
                for (int si = 0; si < 2; ++si) {
                    if (!ctau_freeze[si] && stable_ilam_iters[si] >= 4 &&
                            residualNorm.getValue() > param.rtol.getValue()) {
                        ctau_freeze[si] = true;
                    }
                }
            }

            // NaN-lock detection: sparselinsolve sets dU=0 when the Jacobian has
            // NaN/Inf entries.  A single dU=0 skip is benign; but if the BL state
            // is frozen (no external update can change it), every subsequent
            // iteration produces the same singular Jacobian and the solver spins
            // doing nothing.  Signature: NaN residual, or omega==1.0 with residual
            // not decreasing (dU applied but BL state effectively unchanged).
            {
                double cur = residualNorm.getValue();
                bool stagnant = std::isnan(cur) ||
                    (omega.getValue() == 1.0 && resid_pos > 3 &&
                     prev_resid_val >= 0.0 && cur >= prev_resid_val * 0.9999);
                nan_skip_count = stagnant ? nan_skip_count + 1 : 0;
                prev_resid_val = cur;
                if (nan_skip_count >= 3) {
                    if (failure_mode_out != nullptr)
                        *failure_mode_out = "nan_lock";
                    early_exit = true;
                    break;
                }
            }
        }
    }

    if (!converged && !early_exit && failure_mode_out != nullptr) {
        double rn = (resid_pos > 0) ? resid_buf[(resid_pos - 1) % 8] : 1.0;
        bool had_osc = (stable_ilam_iters[0] >= 15 || stable_ilam_iters[1] >= 15);
        if (had_osc && rn < 1.0)
            *failure_mode_out = "transition_front_oscillation";
        else if (rn >= 1.0)
            *failure_mode_out = "diverged";
        else
            *failure_mode_out = "no_convergence";
    }

    return converged;
}