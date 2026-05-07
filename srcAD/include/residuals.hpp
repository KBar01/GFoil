#pragma once

// get_funcs.hpp must be included before residuals_shared.hpp so get_H etc.
// are visible at instantiation time of the template residual functions.
#include "data_structs.hpp"
#include "get_funcs.hpp"
// Shared physics (upwind, upwind_half, residual_station, residual_station_forced):
#include "residuals_shared.hpp"

template<typename Real> struct Isolc;
template<typename Real> struct Isolv;
template<typename Real> struct Vsol;
template<typename Real> struct Foil;
template<typename Real> struct Param;
template<typename Real> struct Post;
template<typename Real> struct Oper;
template<typename Real> struct Geom;
template<typename Real> struct Wake;
template<typename Real> struct Glob;
template<typename Real> struct Trans;

// AD-specific: wake_sys, residual_transition, residual_transition_forced
// use Vsol<Real>/Foil<Real>/Glob<Real> (template structs from data_structs.hpp)
// which differ from the non-template equivalents in src/include/data_structs.h.

template<typename Real>
void wake_sys(const Vsol<Real>& vsol, const Foil<Real>& foil, const Glob<Real>& glob, const Param<Real>& param,
    Real (&R)[3]) {

    const Real* U = glob.U;

    // --- Get trailing edge states ---
    int il = 0;
    int iu = Ncoords-1;
    int iw = Ncoords;

    Real Ul[4]={0}, Uu[4]={0}, Uw[4]={0};
    for (int i = 0; i < 4; ++i) {
        Ul[i] = U[colMajorIndex(i,il,4)];
        Uu[i] = U[colMajorIndex(i,iu,4)];
        Uw[i] = U[colMajorIndex(i,iw,4)];
    }

    // --- Compute wake shear stress (ctw) ---
    Real ctl, ctu, ctl_Ul[4]={0}, ctu_Uu[4]={0};

    if (vsol.turb[il]) {
        ctl = Ul[2];
        ctl_Ul[0] = 0; ctl_Ul[1] = 0; ctl_Ul[2] = 1; ctl_Ul[3] = 0;
    } else {
        ctl = get_cttr(Ul[0],Ul[1],Ul[2],Ul[3],true,param,ctl_Ul);
    }

    if (vsol.turb[iu]) {
        ctu = Uu[2];
        ctu_Uu[2] = 1;
    } else {
        ctu = get_cttr(Uu[0],Uu[1],Uu[2],Uu[3],true,param,ctu_Uu);
    }

    Real thsum = Ul[0] + Uu[0];
    Real ctw = (ctl * Ul[0] + ctu * Uu[0]) / thsum;

    

    Real ctw_Ul[4]={0}, ctw_Uu[4]={0};

    Real zero = 0;
    for (int i = 0; i < 4; ++i) {
        ctw_Ul[i] = (ctl_Ul[i] * Ul[0] + (i == 0 ? ctl - ctw : zero)) / thsum;
        ctw_Uu[i] = (ctu_Uu[i] * Uu[0] + (i == 0 ? ctu - ctw : zero)) / thsum;
    }

    // --- Residual vector ---
    R[0] = Uw[0] - (Ul[0] + Uu[0]);
    R[1] = Uw[1] - (Ul[1] + Uu[1] + foil.te.hTE);
    R[2] = Uw[2] - ctw;

}

template<typename Real>
void residual_transition(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Real aux1,
    const Real aux2,
    const Param<Real>& param,
    Real (&R)[3]
) {
    const int nNewton = 20;
    const Real dx = x2 - x1;
    
    Real ncrit = param.ncrit;

    // Extract elements of BL States
    Real th1 = U1[0],th2 = U2[0];
    Real ds1 = U1[1],ds2 = U2[1];
    Real sa1 = U1[2],sa2 = U2[2];
    Real ue1 = U1[3],ue2 = U2[3];


    Real damp1, damp1_U1[4];
    
    damp1 = get_damp(th1,ds1,sa1,ue1,param,damp1_U1); // is constant in newton loop
    
    
    Real dampt, dampt_Ut[4]={0};
    Real Rxt, Rxt_xt, dxt;
    Real Ut[4]={0},Ut_xt[4]={0};
    Real w1,w2;
    Real xt = x1 + 0.5 * dx; // initial guess
    
    

    // Not forcing transition, iteratively finding transition location within the panel
    for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
        
        // weights
        w2 = (xt - x1)/dx, w1 =1-w2;

        // State at transition point xt 
        for (int i = 0; i < 4; ++i) {
            Ut[i] = w1*U1[i] + w2*U2[i];
            Ut_xt[i] = (U2[i]-U1[i])/dx;
        }

        Ut[2] = ncrit;   // amplification at transition
        Ut_xt[2] = 0.0;

        dampt = get_damp(Ut[0],Ut[1],Ut[2],Ut[3],param,dampt_Ut);
        dampt_Ut[2] = 0.0;

        Rxt = ncrit - sa1 - 0.5 * (damp1+dampt)*(xt-x1);
        Rxt_xt = -0.5 * (damp1 + dampt);
        for (int i = 0; i < 4; ++i){Rxt_xt -= 0.5*(xt-x1) * dampt_Ut[i]*Ut_xt[i];}

        dxt = -Rxt / Rxt_xt;
        Real dmax = 0.2 * dx * (1.1 - Real(iNewton) / nNewton);
        if (std::abs(dxt) > dmax){ dxt *= dmax / std::abs(dxt);}
        if (std::abs(Rxt) < 1e-10) break;

        xt += dxt;
    }
    
    // Copy to Utl, Utt (laminar and turbulent transition states)
    Real Utl[4] = {0}, Utt[4] = {0};
 
    for (int i = 0; i < 4; ++i) {
        Utl[i] = Ut[i]; Utt[i] = Ut[i]; 
    }
    Utl[2] = ncrit;
    
    Real cttr, cttr_Ut[4]={0};
    cttr = get_cttr(Ut[0],Ut[1],Ut[2],Ut[3],true,param,cttr_Ut);
    Utt[2] = cttr;
    
    //param.turb = false;
    Real Rl[3]={0} ;
    residual_station(U1, Utl, x1, xt, aux1, aux2, false, false, false, param, Rl);
    
    //param.turb = true;
    Real Rt[3]={0}, Rt_U[24]={0}, Rt_x[6]={0};
    residual_station(Utt, U2, xt, x2, aux1, aux2, false, true, false, param, Rt);

    for (int i = 0; i < 3; ++i){
        R[i] = Rl[i] + Rt[i];
    }

}

template<typename Real>
void residual_transition_forced(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Param<Real>& param,
    const Real transPos,
    Real (&R)[3]) {

    const Real dx = x2 - x1;
    Real ncrit = param.ncrit;

    // Extract elements of BL States
    Real th1 = U1[0],th2 = U2[0];
    Real ds1 = U1[1],ds2 = U2[1];
    Real sa1 = U1[2],sa2 = U2[2];
    Real ue1 = U1[3],ue2 = U2[3];
    Real Ut[4]={0}; // state at transiton (initially laminar)
    
    // position is fixed, giving fixed weightings 
    Real xt = transPos;
   
    // Fixed weights and interpolated state
    Real w2 = (xt - x1) / dx;
    Real w1 = 1.0 - w2;

    for (int i = 0; i < 4; ++i){
        Ut[i] = w1*U1[i] + w2*U2[i];
    }
    
    // Copy to Utl, Utt (laminar and turbulent transition states)
    Real Utl[4] = {0}, Utt[4] = {0};

    // coptying lam/turb states at the transition point
    for (int i = 0; i < 4; ++i) {
        Utl[i] = Ut[i];
        Utt[i] = Ut[i];
    }
    
    // initialise the turb state with value of shear stress coeff (note Ut[2] has no effect)
    Real cttr, cttr_Ut[4]={0};
    cttr = get_cttr(Ut[0],Ut[1],Ut[2],Ut[3],true,param,cttr_Ut);
    Utt[2] = cttr;

   
    // residual for laminar section (remember that dont care about amplifcation residual here)
    Real Rl[3]={0};
    residual_station(U1, Utl, x1, xt, Real(0.0), Real(0.0), false, false, false, param, Rl);
    Rl[2] = 0.0;
    
    // residual for the turbulent section, here residual for shear stress matters
    Real Rt[3]={0};
    residual_station(Utt, U2, xt, x2, Real(0.0), Real(0.0), false, true, false, param, Rt);

    // sum residuals for total panel residuals
    for (int i = 0; i < 3; ++i){
        R[i] = Rl[i] + Rt[i];
    }

}
