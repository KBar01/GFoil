#pragma once

#include <cmath>
// Callers must include real_type.h / real_type.hpp, and the shared
// get_funcs.hpp (for get_cf, get_uk) before including this header.
#include "get_funcs.hpp"




// Performs cubic interpolation at x using 4 surrounding points (xs, ys)
template<typename Real>
Real cubic_interp(Real x, const Real* xs, const Real* ys) {
    Real result = 0.0;
    for (int i = 0; i < 4; ++i) {
        Real term = ys[i];
        for (int j = 0; j < 4; ++j) {
            if (i != j)
                term *= (x - xs[j]) / (xs[i] - xs[j]);
        }
        result += term;
    }
    return result;
}

template<typename Real>
Real linear_interp(Real x, const Real* xs, const Real* ys) {
    return ys[0] + (ys[1] - ys[0]) * (x - xs[0]) / (xs[1] - xs[0]);
}

template<typename Real>
Real quadratic_interp(Real x, const Real* xs, const Real* ys) {
    Real result = 0.0;
    for (int i = 0; i < 3; ++i) {
        Real term = ys[i];
        for (int j = 0; j < 3; ++j) {
            if (i != j)
                term *= (x - xs[j]) / (xs[i] - xs[j]);
        }
        result += term;
    }
    return result;
}

template<typename Real>
Real adaptive_interp(Real x, const Real* xs, const Real* ys, int n) {
    if (n == 4)
        return cubic_interp(x, xs, ys);
    else if (n == 3)
        return quadratic_interp(x, xs, ys);
    else if (n == 2)
        return linear_interp(x, xs, ys);
    else if (n == 1)
        return ys[0];
    else
        return Real(0.0); // n==0: laminar surface; caller checks tauMax before use
}

template<typename Real>
int find_interp_position(const Real* xcoords, int start, int end, Real x_target) {

    // For top surface : finds node just before sampling pos
    // For bot surface : finds node just before sampling pos
    
    int found_idx = -1;
    bool is_increasing = xcoords[start] < xcoords[start + 1];

    for (int i = start; i < end - 1; ++i) {
        if (is_increasing) {
            if (xcoords[i] <= x_target && x_target < xcoords[i + 1]) {
                found_idx = i;
                break;
            }
        } else {
            if (xcoords[i] >= x_target && x_target > xcoords[i + 1]) {
                found_idx = i+1;
                break;
            }
        }
    }

    return found_idx;
}

template<typename Real, typename TurbT>
void get_nodes(int topFoundIdx, int botFoundIdx, Real x_target, int* topNodeList, int &topNnodes, int* botNodeList, int &botNnodes, const TurbT* turb) {

    // Top surface (increasing x, LE→TE): build stencil of up to 4 turbulent nodes.
    // topFoundIdx is the node just before x_target; start one node earlier for symmetry.
    int topStart = topFoundIdx - 1;
    if (topFoundIdx < 0 || !turb[topFoundIdx]) {
        topNnodes = 0;
    } else {
        if (topStart < 0 || !turb[topStart]) { topStart = topFoundIdx; }
        topNnodes = (Ncoords - 1 - topStart) + 1;
        if (topNnodes > 4) { topNnodes = 4; }
    }
    for (int i = 0; i < topNnodes; ++i) { topNodeList[i] = topStart + i; }

    // Bot surface (decreasing x, TE→LE): build stencil of up to 4 turbulent nodes.
    // botFoundIdx is the node just before x_target; start one node further toward LE.
    int botStart = botFoundIdx + 1;
    if (botFoundIdx < 0 || !turb[botFoundIdx]) {
        botNnodes = 0;
    } else {
        if (!turb[botStart]) { botStart = botFoundIdx; }
        botNnodes = botStart + 1;  // nodes botStart..0 = botStart+1 available
        if (botNnodes > 4) { botNnodes = 4; }
    }
    for (int i = 0; i < botNnodes; ++i) { botNodeList[i] = botStart - i; }
}

template<typename Real>
void interp_BL_states(const int* topIdx,const int* botIdx, const int topNnodes, const int botNnodes, const Real x_target, const Real* xcoords, const Real* states, Real* topInterpStates, Real* botInterpStates){



    // Do top surface first 
    Real txs[4] = {
        xcoords[topIdx[0]],
        xcoords[topIdx[1]],
        xcoords[topIdx[2]],
        xcoords[topIdx[3]]
    };

    for (int q = 0; q < 4; ++q) {
        Real tys[4] = {
            states[colMajorIndex(q,topIdx[0],4)],
            states[colMajorIndex(q,topIdx[1],4)],
            states[colMajorIndex(q,topIdx[2],4)],
            states[colMajorIndex(q,topIdx[3],4)]
        };
        topInterpStates[q] = adaptive_interp(x_target,txs,tys,topNnodes);
    }

    // then bottom surface

    Real bxs[4] = {
        xcoords[botIdx[0]],
        xcoords[botIdx[1]],
        xcoords[botIdx[2]],
        xcoords[botIdx[3]]
    };

    for (int q = 0; q < 4; ++q) {
        Real bys[4] = {
            states[colMajorIndex(q,botIdx[0],4)],
            states[colMajorIndex(q,botIdx[1],4)],
            states[colMajorIndex(q,botIdx[2],4)],
            states[colMajorIndex(q,botIdx[3],4)]
        };
        botInterpStates[q] = adaptive_interp(x_target,bxs,bys,botNnodes);
    }
}


// Computes dp/dx at x_target using cubic interpolation + central finite difference
template<typename Real, typename OperT>
Real interpolate_dpdx(const Real* xcoords, const Real* Cps, const int* nodeIdx, const int nodeN, Real x_target, const OperT& oper,const Real chordScale,const Real Uinf) {
    
    
    Real h = 1e-6;  // unit-chord step; differentiates the cubic analytically to machine precision

    // Get x positions slightly left and right of target
    Real x_plus  = x_target + h;
    Real x_minus = x_target - h;

    Real xs[4] = {
        xcoords[nodeIdx[0]],
        xcoords[nodeIdx[1]],
        xcoords[nodeIdx[2]],
        xcoords[nodeIdx[3]]
    };

    Real ps[4] = {
        Cps[nodeIdx[0]],
        Cps[nodeIdx[1]],
        Cps[nodeIdx[2]],
        Cps[nodeIdx[3]]
    };

    // Interpolate pressure at x+h and x-h
    Real CpPlus  = adaptive_interp(x_plus,xs,ps,nodeN);
    Real CpMinus = adaptive_interp(x_minus,xs,ps,nodeN);

    // Central finite difference
    Real dpdx = ((CpPlus - CpMinus) / (2.0*(h*chordScale))) * (0.5 * oper.rho * Uinf*Uinf);

    return dpdx;
}

template<typename Real, typename TurbT, typename ParamT>
Real interpolate_cf(const Real* xcoords, const Real* states, const int* nodeIdx, const int nodeN, Real x_target, const TurbT* turb, const ParamT& param) {
    
    Real xs[4] = {
        xcoords[nodeIdx[0]],
        xcoords[nodeIdx[1]],
        xcoords[nodeIdx[2]],
        xcoords[nodeIdx[3]]
    };


    Real cf[4];

    int indexes[4] = {nodeIdx[0], nodeIdx[1],nodeIdx[2], nodeIdx[3]} ;
    Real cf_U[4] = {0};
    for (int i=0;i<4;++i){
        cf[i] = get_cf(
            states[colMajorIndex(0,indexes[i],4)],
            states[colMajorIndex(1,indexes[i],4)],
            states[colMajorIndex(2,indexes[i],4)],
            states[colMajorIndex(3,indexes[i],4)],
            turb[indexes[i]],
            false,
            param,
            cf_U
        );
    }
    // Interpolate 
    Real Cf95  = adaptive_interp(x_target,xs,cf,nodeN);

    return Cf95;
}

template<typename Real, typename OperT, typename TurbT, typename ParamT>
void interpolate_at_95_both_surfaces(const Real* xcoords, const Real* states, const Real* Cps, const OperT& oper, const TurbT* turb, const ParamT& param,
    Real (&topBLStates)[7], Real (&botBLStates)[7], const Real Uinf, const Real x_target, const Real chordScale) {

    /* State order: theta, delta*, tau_max, Ue, dpdx, tau_wall, delta99 */

    int foundIndexBot = find_interp_position(xcoords, 0, 98, x_target);
    int foundIndexTop = find_interp_position(xcoords, Ncoords - 98, Ncoords - 1, x_target);

    int topIdx[4] = {0}, botIdx[4] = {0}, topN, botN;
    get_nodes(foundIndexTop, foundIndexBot, x_target, topIdx, topN, botIdx, botN, turb);

    interp_BL_states(topIdx, botIdx, topN, botN, x_target, xcoords, states, topBLStates, botBLStates);

    Real dpdxBot = interpolate_dpdx(xcoords, Cps, botIdx, botN, x_target, oper, chordScale, Uinf);
    Real dpdxTop = interpolate_dpdx(xcoords, Cps, topIdx, topN, x_target, oper, chordScale, Uinf);

    // ── Bottom surface: convert to dimensional quantities ────────────────────
    Real ignore;
    Real UeCorrected  = get_uk(botBLStates[3], param, ignore) * Uinf;
    botBLStates[3]    = UeCorrected;
    botBLStates[2]    = (botBLStates[2] * botBLStates[2]) * (oper.rho * botBLStates[3] * botBLStates[3]);
    botBLStates[4]    = dpdxBot;
    Real cfBot        = interpolate_cf(xcoords, states, botIdx, botN, x_target, turb, param);
    Real tauWallBot   = (cfBot / 2) * oper.rho * botBLStates[3] * botBLStates[3];
    botBLStates[5]    = tauWallBot;
    botBLStates[0]   *= chordScale;
    botBLStates[1]   *= chordScale;
    {
        Real H_bot = (botBLStates[0] > Real(0.0)) ? botBLStates[1] / botBLStates[0] : Real(2.0);
        botBLStates[6] = (tauWallBot != 0.0 && H_bot > Real(1.0) + Real(1e-4))
            ? botBLStates[0] * (3.15 + 1.72 / (H_bot - Real(1.0))) + botBLStates[1]
            : Real(0.0);
    }

    // ── Top surface: convert to dimensional quantities ───────────────────────
    UeCorrected       = get_uk(topBLStates[3], param, ignore) * Uinf;
    topBLStates[3]    = UeCorrected;
    topBLStates[2]    = (topBLStates[2] * topBLStates[2]) * (oper.rho * topBLStates[3] * topBLStates[3]);
    topBLStates[4]    = dpdxTop;
    Real cfTop        = interpolate_cf(xcoords, states, topIdx, topN, x_target, turb, param);
    Real tauWallTop   = (cfTop / 2) * oper.rho * topBLStates[3] * topBLStates[3];
    topBLStates[5]    = tauWallTop;
    topBLStates[0]   *= chordScale;
    topBLStates[1]   *= chordScale;
    {
        Real H_top = (topBLStates[0] > Real(0.0)) ? topBLStates[1] / topBLStates[0] : Real(2.0);
        topBLStates[6] = (tauWallTop != 0.0 && H_top > Real(1.0) + Real(1e-4))
            ? topBLStates[0] * (3.15 + 1.72 / (H_top - Real(1.0))) + topBLStates[1]
            : Real(0.0);
    }
}



