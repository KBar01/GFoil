#pragma once

// extract_BL_TE.hpp — sample the trailing-edge BL state for the acoustic model.
//
// interpolate_at_95_both_surfaces() interpolates the converged boundary-layer
// quantities to the acoustic sampling station near the trailing edge and packs
// the two 7-slot vectors [theta, deltaStar, tauMax, Ue, dpdx, tauWall, delta]
// (upper/lower) that feed calc_OASPL -> calc_WPS. Bridges the aero solve and the
// acoustics (sound.hpp). Symbols: see NOMENCLATURE.md.
#include <cmath>
#include <stdexcept>
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
    else
        return 0.0; // fallback error
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
void get_nodes(int topFoundIdx,int botFoundIdx, Real x_target, int* topNodeList,int &topNnodes, int* botNodeList,int &botNnodes, const TurbT* turb){

    // do top surface :
    int topStart = topFoundIdx - 1 ; // Ideal starting position for cubic interp to have 2 nodes either side of sampling position
    topNnodes = 4 ;
    if (turb[topFoundIdx] == false){
        topNnodes = 0;
    }
    else{
        if (turb[topStart] == false){topStart += 1; } // original start is not turbulent, shift start down a node
        topNnodes = (Ncoords-1 - topStart) + 1 ;  // how many available nodes to use for interp
        if (topNnodes>4){topNnodes=4;} // limit to 4 nodes for cubic
    }

    for (int i=0;i<topNnodes;++i){topNodeList[i] = topStart+i ;}


    // do bot surface :
    int botStart = botFoundIdx + 1 ; // Ideal starting position for cubic interp to have 2 nodes either side of sampling position
    botNnodes = 4;
    if (turb[botFoundIdx] == false){
        botNnodes = 0;
    }
    else{
        if (turb[botStart] == false){botStart -= 1; } // original start is not turbulent, shift start down a node
        botNnodes = botStart ;  // how many available nodes to use for interp (write the OUT param, not a shadow)
        if (botNnodes>4){botNnodes=4;} // limit to 4 nodes for cubic
    }

    for (int i=0;i<botNnodes;++i){botNodeList[i] = botStart-i ;}
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
    
    
    Real h = 1e-6 ;  // Small step size for derivative approximation TODO: verfiy step is correct (convergence)

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

// Single-point TE BL sampler. Fills the two fully post-processed 7-slot vectors
// [theta, delta*, tau_max, Ue, dpdx, tau_wall, delta99] at chord fraction x_target.
// The public interpolate_at_95_both_surfaces() dispatcher below calls it (once for
// a scalar TE sample, N times for a window). x_target must be < 1.0: the trailing-
// edge node is degenerate for the interpolation stencil (find_interp_position has
// no bracket to return), which is enforced by the dispatcher. (The former x_target
// == 1.0 special case averaged over a hardcoded 0.96–0.985 window — removed now
// that an explicit [x_lo,x_hi] window can be requested directly.)
template<typename Real, typename OperT, typename TurbT, typename ParamT>
void interpolate_BL_single(const Real* xcoords, const Real* states, const Real*Cps, const OperT& oper, const TurbT* turb, const ParamT& param,
    Real (&topBLStates)[7],Real (&botBLStates)[7],const Real Uinf, const Real x_target, const Real chordScale) {


    /* State order: theta, delta*, tau_max, Ue, dpdx, tau_wall, delta 99% thickness*/
    
    // find index of node before sampling position (top and bottom)
    int foundIndexBot = find_interp_position(xcoords,0,98,x_target);
    int foundIndexTop = find_interp_position(xcoords,Ncoords-98, Ncoords-1,x_target);

    // find the indexes to to the interpolation over 
    int topIdx[4] = {0},botIdx[4] = {0}, topN, botN ;
    get_nodes(foundIndexTop,foundIndexBot,x_target,topIdx,topN,botIdx,botN,turb);

    interp_BL_states(topIdx,botIdx,topN,botN,x_target,xcoords,states,topBLStates,botBLStates);
    
    Real dpdxBot = interpolate_dpdx(xcoords,Cps,botIdx,botN,x_target,oper,chordScale,Uinf);
    Real dpdxTop = interpolate_dpdx(xcoords,Cps,topIdx,topN,x_target,oper,chordScale,Uinf);
    
    // ---------------------- bottom surface dimensionals -----------------------------------------
    Real ignore;
    Real UeCorrected = (get_uk(botBLStates[3],param,ignore)) * Uinf;
    botBLStates[3] = UeCorrected;
    // BLstate is C_tau ^ 0.5 , and C_tau = tau_max / (rho * Ue^2)
    Real tauMaxBot = (botBLStates[2] * botBLStates[2]) * (oper.rho * (botBLStates[3]*botBLStates[3])) ;
    
    botBLStates[2] = tauMaxBot;
    botBLStates[4] = dpdxBot;

    // tau_wall = Cf*(rho * ue^2) / 2
    Real cfBot = interpolate_cf(xcoords,states,botIdx,botN,x_target,turb,param);
    Real tauWallBot = (cfBot/2) * oper.rho * botBLStates[3] * botBLStates[3] ;
    botBLStates[5] = tauWallBot ;
    // scaling theta and delta* by given chord 
    botBLStates[0] *= chordScale ;
    botBLStates[1] *= chordScale ;


    // now get 99% thickness
    Real deltaBot = 0.0;
    Real frictionVel = std::sqrt(tauWallBot/oper.rho) ;
    if (tauWallBot!=0.0){
        deltaBot = botBLStates[0]*(3.15 + 1.72/((botBLStates[1]/botBLStates[0]) - 1)) + botBLStates[1] ;
    }
    botBLStates[6] = deltaBot;

    // --------------------------------- top surface dimensionals ----------------------------------------
    UeCorrected = (get_uk(topBLStates[3],param,ignore)) * Uinf;
    topBLStates[3] = UeCorrected;
    
    // BLstate is C_tau ^ 0.5 , and C_tau = tau_max / (rho * Ue^2)
    Real tauMaxTop = (topBLStates[2] * topBLStates[2]) * (oper.rho * (topBLStates[3]*topBLStates[3])) ;

    topBLStates[2] = tauMaxTop;
    topBLStates[4] = dpdxTop;
    
    // tau_wall = Cf*(rho * ue^2) / 2
    Real cfTop = interpolate_cf(xcoords,states,topIdx,topN,x_target,turb,param);
    Real tauWallTop = (cfTop/2) * oper.rho * topBLStates[3] * topBLStates[3] ;
    topBLStates[5] = tauWallTop ;


    topBLStates[0] *= chordScale ;
    topBLStates[1] *= chordScale ;
    // now get 99% thickness
    Real deltaTop = 0.0;
    frictionVel = std::sqrt(tauWallTop/oper.rho) ;
    if (tauWallTop!=0.0){
        deltaTop = topBLStates[0]*(3.15 + 1.72/((topBLStates[1]/topBLStates[0]) - 1)) + topBLStates[1] ;
    }
    
    topBLStates[6] = deltaTop;

}


// Number of trapezoidal quadrature stations spanning a BL-averaging window.
// 9 gives good trapezoidal accuracy at modest extra tape cost; the endpoints are
// included so x_lo / x_hi are sampled (interpolated) directly.
#ifndef NWINDOW_SAMPLES
#define NWINDOW_SAMPLES 9
#endif

// Backward-compatible TE sampling entry point. All sample points must be < 1.0.
//   x_hi <= x_lo : single-point sample at x_lo (byte-identical to the legacy path).
//   x_hi  > x_lo : Option-A BL-averaged window. The fully post-processed 7-slot
//                  BL/WPS input vectors are evaluated at NWINDOW_SAMPLES uniformly
//                  spaced stations across [x_lo, x_hi] (endpoints interpolated, not
//                  snapped) and trapezoidally averaged in x/c before a single
//                  downstream Amiet evaluation. This removes the single-node
//                  wall-pressure lever the optimiser was exploiting: a smooth
//                  one-node surface undulation can no longer swing the predicted
//                  OASPL because the WPS input is integrated over the window.
//
// AD note: x_lo/x_hi are passive (cast from a double input, never registered), so
// the station positions xs and trapezoid weights w are passive constants. The
// sampled BL quantities (tmp*) are taped through glob.U and the node x-positions,
// so the average remains correctly differentiated w.r.t. y / alpha. All-Real
// arithmetic — no fabs/hypot/getValue here.
template<typename Real, typename OperT, typename TurbT, typename ParamT>
void interpolate_at_95_both_surfaces(const Real* xcoords, const Real* states, const Real* Cps,
    const OperT& oper, const TurbT* turb, const ParamT& param,
    Real (&topBLStates)[7], Real (&botBLStates)[7], const Real Uinf,
    const Real x_lo, const Real x_hi, const Real chordScale) {

    // No sample point may reach the trailing-edge node (x/c == 1.0): it is
    // degenerate for the interpolation stencil (find_interp_position returns no
    // bracket). The legacy x_target==1.0 averaging hack that papered over this was
    // removed — request an explicit window strictly inside (0,1) for TE-region
    // sampling. The scalar point is x_lo; the top window station is x_hi. (The
    // Python layer also enforces this; the C++ guard covers the standalone/JSON
    // path.)
    if (x_lo >= 1.0 || x_hi >= 1.0) {
        throw std::invalid_argument(
            "TEsample x/c must be < 1.0 (the trailing-edge node is degenerate for "
            "interpolation; the legacy 1.0 averaging was removed — use a window "
            "[x_lo,x_hi] inside (0,1) for trailing-edge sampling)");
    }

    // Scalar / degenerate-window path: reproduce the legacy single-point result exactly.
    if (!(x_hi > x_lo)) {
        interpolate_BL_single(xcoords, states, Cps, oper, turb, param,
                              topBLStates, botBLStates, Uinf, x_lo, chordScale);
        return;
    }

    // Windowed Option-A path: trapezoidal average in x/c of the single-point inputs.
    // Average = (1/width) * trapz(f) ; the window width cancels, leaving station
    // weights [0.5, 1, ..., 1, 0.5] / (N-1).
    constexpr int N = NWINDOW_SAMPLES;
    Real topAccum[7] = {0.0}, botAccum[7] = {0.0};
    for (int i = 0; i < N; ++i) {
        Real frac = static_cast<Real>(i) / static_cast<Real>(N - 1);
        Real xs   = x_lo + (x_hi - x_lo) * frac;
        Real w    = ((i == 0 || i == N - 1) ? 0.5 : 1.0) / static_cast<Real>(N - 1);

        Real tmpTop[7], tmpBot[7];
        interpolate_BL_single(xcoords, states, Cps, oper, turb, param,
                              tmpTop, tmpBot, Uinf, xs, chordScale);
        for (int k = 0; k < 7; ++k) {
            topAccum[k] += w * tmpTop[k];
            botAccum[k] += w * tmpBot[k];
        }
    }
    for (int k = 0; k < 7; ++k) {
        topBLStates[k] = topAccum[k];
        botBLStates[k] = botAccum[k];
    }
}



