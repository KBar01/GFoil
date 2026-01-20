#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"


void update_state(const Oper&oper, const Param&param, Glob&glob,Vsol&vsol) {
    
    constexpr int Nsys = Ncoords+Nwake;

    

    // Check for imaginary components (not directly applicable in C++, skip unless you handle complex numbers)
    // max ctau
    std::vector<int> It;
    Real ctmax = -1e20;
    for (int i = 0; i < Nsys; ++i) {
        if (vsol.turb[i]) {
            It.push_back(i);
            Real ct = glob.U[colMajorIndex(2, i, 4)];
            if (ct > ctmax){ ctmax = ct;}
        }
    }

    Real omega = 1.0;

    // limit theta and delta*
    const Real one = 1;
    for (int k = 0; k < 2; ++k) {
        Real fmin = 1e20;
        for (int i = 0; i < Nsys; ++i) {
            Real Uk = glob.U[colMajorIndex(k, i, 4)];
            Real dUk = glob.dU[colMajorIndex(k, i, 4)];
            if (Uk != 0.0) {
                Real ratio = dUk / Uk;
                if (ratio < fmin) fmin = ratio;
            }
        }
        Real om = (fmin < -0.5) ? std::abs(0.5 / fmin) : one;
        if (om < omega) {omega = om;}
    }

    // limit negative amp/ctau
    for (int i = 0; i < Nsys; ++i) {
        Real Uk = glob.U[colMajorIndex(2, i, 4)];
        Real dUk = glob.dU[colMajorIndex(2, i, 4)];
        if (!vsol.turb[i] && Uk < 0.2) continue;
        if (vsol.turb[i] && Uk < 0.1 * ctmax) continue;
        if (Uk == 0.0 || dUk == 0.0) continue;
        if ((Uk + dUk) < 0.0) {
            Real om = 0.8 * std::abs(Uk / dUk);
            if (om < omega) {omega = om;}
        }
    }

    // prevent big amp changes (non-turbulent)
    Real dumax = 0.0;
    for (int i = 0; i < Nsys; ++i) {
        if (!vsol.turb[i]) {
            Real dUk = glob.dU[colMajorIndex(2, i, 4)];
            if (std::abs(dUk) > dumax) dumax = std::abs(dUk);
        }
    }
    Real om  = (dumax > 0) ? std::abs(2.0 / dumax) : one;
    if (om < omega) {omega = om;}


    // prevent big ctau changes (turbulent)
    dumax = 0.0;
    for (int i = 0; i < Nsys; ++i) {
        if (vsol.turb[i]) {
            Real dUk = glob.dU[colMajorIndex(2, i, 4)];
            if (std::abs(dUk) > dumax) dumax = std::abs(dUk);
        }
    }
    om  = (dumax > 0) ? std::abs(0.05 / dumax) : one;
    if (om < omega) {omega = om;}
    

    // prevent large ue changes
    dumax = 0.0;
    for (int i = 0; i < Nsys; ++i) {
        Real dUk = glob.dU[colMajorIndex(3, i, 4)];
        Real ratio = std::abs(dUk / oper.Vinf);
        if (ratio > dumax) dumax = ratio;
    }
    om  = (dumax > 0) ? 0.2/dumax : one;
    if (om < omega) {omega = om;}
    

    // take update
    for (int i = 0; i < 4*Nsys; ++i) {
        glob.U[i] += omega * glob.dU[i];
    }

    // fix bad Hk
    for (int si = 0; si < 3; ++si) {
        
        Real Hkmin = (si == 2) ? 1.00005 : 1.02;
        
        const std::vector<int>& Is = vsol.Is[si];
        const int nSurfPoints = Is.size();

        for (int ii = 0; ii < nSurfPoints; ++ii) {
            
            int j = Is[ii];
            Real Uj[4] = {
                glob.U[colMajorIndex(0, j, 4)],
                glob.U[colMajorIndex(1, j, 4)],
                glob.U[colMajorIndex(2, j, 4)],
                glob.U[colMajorIndex(3, j, 4)]
            };

            Real Hk, Hk_U[4];
            Hk = get_Hk(Uj[0],Uj[1],Uj[3],param,Hk_U);
            if (Hk < Hkmin) {
                glob.U[colMajorIndex(1, j, 4)] += 2.0 * (Hkmin - Hk) * glob.U[colMajorIndex(1, j, 4)];
            }
        }
    }

    // fix negative ctau
    for (int i : It) {
        if (glob.U[colMajorIndex(2, i, 4)] < 0.0) {
            glob.U[colMajorIndex(2, i, 4)] = 0.1 * ctmax;
        }
    }
}
