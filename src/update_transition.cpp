#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.hpp"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.hpp"
#include "vector_ops.hpp"

int march_amplification(Glob &glob, Vsol &vsol, Isol &isol, int si, const Param&param,
                        Real* amp_break = nullptr, bool* forced_break = nullptr) {

    const std::vector<int> &Is = vsol.Is[si];
    int N = Is.size();

    glob.U[colMajorIndex(2,Is[0],4)] = 0.0; // initial amplification

    int i = 1;
    while (i < N) {

        int i1 = Is[i-1];
        int i2 = Is[i];

        Real U1[4], U2[4];
        for (int j = 0; j < 4; ++j) {
            U1[j] = glob.U[colMajorIndex(j,i1,4)];
            U2[j] = glob.U[colMajorIndex(j,i2,4)];
        }
        if (vsol.turb[i2]) {U2[2] = U1[2]*1.01;}

        Real dx = isol.distFromStag[i2] - isol.distFromStag[i1];

        constexpr int nNewton = 20;
        const Real one = 1, zero=0;

        for (int iNewton=0; iNewton < nNewton; ++iNewton) {

            Real damp1, damp2, damp1_U[4], damp2_U[4];
            damp1 = get_damp(U1[0],U1[1],U1[2],U1[3],param,damp1_U);
            damp2 = get_damp(U2[0],U2[1],U2[2],U2[3],param,damp2_U);

            Real damp, damp_U[8];
            damp = upwind_half(damp1, damp1_U, damp2, damp2_U, damp_U);

            Real Ramp = U2[2] - U1[2] - damp*dx;
            if (std::abs(Ramp) < 1e-12) break;

            Real Ramp_U[8] = {0.0};
            Ramp_U[2] = -1.0;
            Ramp_U[6] = 1.0;
            for (int j = 0; j < 8; ++j){ Ramp_U[j] -= damp_U[j]*dx;}

            Real dU = -Ramp / Ramp_U[6];
            Real dmax = 0.5 * (1.01 - static_cast<Real>(iNewton)/nNewton);
            Real omega = (std::abs(dU) > dmax) ? dmax/std::abs(dU) : one;
            U2[2] += omega*dU;
        }

        // Check 1 — natural transition (ncrit wins, regardless of forced)
        if (U2[2] > param.ncrit) {
            if (amp_break) *amp_break = U2[2];
            break;
        }

        // Check 2 — forced transition: xift falls in this interval
        if (param.forcet && param.xift > 0.0) {
            double xi1_val = isol.distFromStag[i1].getValue();
            double xi2_val = isol.distFromStag[i2].getValue();
            if (xi1_val <= param.xift && param.xift < xi2_val) {
                if (amp_break)    *amp_break    = U2[2];
                if (forced_break) *forced_break = true;
                break;
            }
        }

        glob.U[colMajorIndex(2,i2,4)] = U2[2];
        ++i;
    }

    return i - 1;
}


void update_transition(Glob &glob, Vsol &vsol, Isol &isol, Param &param,
                       const Foil& foil,
                       int newtonIter,
                       bool freeze_ctau_top, bool freeze_ctau_bot,
                       double* prev_amp_top,  double* prev_amp_bot,
                       double* prev_ctau_top, double* prev_ctau_bot) {

    for (int si = 0; si < 2; ++si) {

        const std::vector<int> &Is = vsol.Is[si];
        int nSurfPoints = Is.size();

        // Precompute forced transition arc-length for this surface.
        // update_transition is forward-only; .getValue() is safe here.
        vsol.forcet[si] = false;
        vsol.xift[si]   = 0.0;
        if (param.xft_xc[si] < 1.0 - 1e-9) {
            double x_max = 0.0;
            for (int k = 0; k < Ncoords; ++k)
                x_max = std::max(x_max, foil.x[2*k].getValue());
            double xft_abs = param.xft_xc[si] * x_max;

            for (int k = 1; k < nSurfPoints; ++k) {
                double x_prev = foil.x[2 * Is[k-1]].getValue();
                double x_curr = foil.x[2 * Is[k  ]].getValue();
                if ((x_prev - xft_abs) * (x_curr - xft_abs) <= 0.0) {
                    double xi_prev = isol.distFromStag[Is[k-1]].getValue();
                    double xi_curr = isol.distFromStag[Is[k  ]].getValue();
                    double frac = (x_curr == x_prev) ? 0.0 :
                                  (xft_abs - x_prev) / (x_curr - x_prev);
                    vsol.xift[si]   = xi_prev + (xi_curr - xi_prev) * frac;
                    vsol.forcet[si] = true;
                    break;
                }
            }
        }
        param.forcet = vsol.forcet[si];
        param.xift   = vsol.xift[si];

        // find current last laminar station
        int ilam0 = nSurfPoints - 1;
        for (int i = 0; i < nSurfPoints; ++i) {
            if (vsol.turb[Is[i]]) {
                ilam0 = i - 1;
                break;
            }
        }

        // copy current amp/shear
        Real sa[Ncoords];
        for (int state=0;state<Ncoords;++state){
            sa[state] = glob.U[colMajorIndex(2,state,4)];
        }

        Real amp_break = 0.0;
        bool was_forced_break = false;
        int ilam = march_amplification(glob, vsol, isol, si, param, &amp_break, &was_forced_break);

        // Apply jump cap for advance direction on both natural and forced transition.
        // The cap prevents large ctau-state discontinuities that cause NaN in the Jacobian.
        // For forced transition, use a larger cap so the target is reached quickly.
        // For forced transition, hysteresis is skipped (location is geometrically fixed).
        if (ilam < ilam0) {
            int max_advance;
            if (was_forced_break) {
                // Reach the forced target in ~4 steps regardless of newtonIter
                max_advance = std::max(1, (ilam0 - ilam + 3) / 4);
            } else {
                max_advance = (newtonIter < 5) ? 1 : 3;
            }
            ilam = std::max(ilam, ilam0 - max_advance);
        } else if (!was_forced_break && ilam > ilam0 &&
                   (ilam - ilam0 == 1) && (ilam0 + 1 < nSurfPoints)) {
            // Hysteresis: suppress spurious 1-node retreat for free transition only.
            Real amp_first_turb = glob.U[colMajorIndex(2, Is[ilam0+1], 4)];
            if (amp_first_turb >= param.ncrit - param.ncrithyst) {
                ilam = ilam0;
            }
        }

        if (ilam == ilam0) {
            // Keep march-computed laminar amps (they satisfy the eN ODE exactly).
            // Only restore turbulent nodes whose ctau march overwrote with amp values.
            for (int state = 0; state < Ncoords; ++state) {
                if (vsol.turb[state]) {
                    glob.U[colMajorIndex(2, state, 4)] = sa[state];
                }
            }

            // ctau freeze: when a period-N limit cycle is detected by solve_coupled,
            // anchor the first turbulent node's ctau to equilibrium (get_cttr) instead
            // of the Newton-updated value that is cycling.  This is additive — it runs
            // after the restore loop above and overrides only the single transition-front
            // node.  Only fires when ilam0+1 is a valid turbulent node.
            bool freeze = (si == 0) ? freeze_ctau_bot : freeze_ctau_top;
            if (freeze) {
                // Average ctau at the first turbulent node (Is[ilam0+1]) between
                // the current Newton-updated value and the previous Newton-updated
                // value.  Storing the pre-averaging Newton value (not the averaged
                // result) means a period-2 cycle (alternating C_A / C_B) collapses
                // to (C_A+C_B)/2 after just two freeze iterations.
                // Note: get_cttr was tried here but it reads the oscillating BL state
                // and therefore itself oscillates; plain averaging is more robust.
                double* prev_ctau = (si == 0) ? prev_ctau_bot : prev_ctau_top;
                if (prev_ctau != nullptr && (ilam0 + 1 < nSurfPoints) &&
                        vsol.turb[Is[ilam0 + 1]]) {
                    double cn = glob.U[colMajorIndex(2, Is[ilam0 + 1], 4)].getValue();
                    if (prev_ctau[0] >= 0.0) {
                        glob.U[colMajorIndex(2, Is[ilam0 + 1], 4)] =
                            Real((prev_ctau[0] + cn) * 0.5);
                    }
                    prev_ctau[0] = cn;
                }
            }

            continue;
        }

        if (ilam < ilam0) {
            bool turb = true;
            Real sa0, cttr_U[4];
            sa0 = get_cttr(glob.U[colMajorIndex(0,Is[ilam+1],4)],
                glob.U[colMajorIndex(1,Is[ilam+1],4)],
                glob.U[colMajorIndex(2,Is[ilam+1],4)],
                glob.U[colMajorIndex(3,Is[ilam+1],4)],
                turb,param,cttr_U);

            Real sa1 = (ilam0 < nSurfPoints-1) ? glob.U[colMajorIndex(2,Is[ilam0+1],4)] : sa0;

            const Real zero=0,one=1;
            Real xi_start = isol.distFromStag[Is[ilam+1]];
            Real xi_end = isol.distFromStag[Is[std::min(ilam0+1, nSurfPoints-1)]];
            Real dx = xi_end - xi_start;

            for (int i = ilam+1; i <= ilam0; ++i) {
                Real f = (dx == 0 || i == ilam+1) ? zero : (isol.distFromStag[Is[i]]-xi_start)/dx;
                if ((ilam+1) == ilam0) f = one;
                Real sa_interp = sa0 + f*(sa1 - sa0);
                glob.U[colMajorIndex(2,Is[i],4)] = sa_interp;
                vsol.turb[Is[i]] = true;
            }
        }
        else if (ilam > ilam0){
            for (int i = ilam0; i <= ilam; ++i)
                vsol.turb[Is[i]] = false;
        }
    }
}
