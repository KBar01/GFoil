#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"
#include "residuals.h"  // pulls in residuals.hpp + non-template wrappers for internal calls


void wake_sys(const Vsol& vsol, const Foil& foil, const Glob& glob, const Param& param,
    Real (&R)[3], Real (&R_U)[36], int (&J)[3]) {

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

    
    // --- Jacobian block indices ---
    J[0] = il;
    J[1] = iu;
    J[2] = iw;

    // --- Residual Jacobian matrix: 3x12 ---
    // R_Ul (cols 0-3): -I (2 rows) and -ctw_Ul
    R_U[0] = -1.0;
    R_U[colMajorIndex(1,1,3)] = -1.0;
    for (int i = 0; i < 4; ++i){
        R_U[colMajorIndex(2,i,3)] = -ctw_Ul[i];
    }

    // R_Uu (cols 4-7): -I (2 rows) and -ctw_Uu
    R_U[colMajorIndex(0,4,3)] = -1.0;
    R_U[colMajorIndex(1,5,3)] = -1.0;
    for (int i = 0; i < 4; ++i){
        R_U[colMajorIndex(2,i+4,3)] = -ctw_Uu[i];
    }
    // R_Uw (cols 8-11): identity
    R_U[colMajorIndex(0,8,3)]  = 1.0;
    R_U[colMajorIndex(1,9,3)]  = 1.0;
    R_U[colMajorIndex(2,10,3)] = 1.0;
}

void residual_transition(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Real aux1,
    const Real aux2,
    const Param& param,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6]
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
    
    
    Real Rxt_U[8] = {0};
    for (int i = 0; i < 4; ++i) {
        Rxt_U[i]   = -0.5 * (xt-x1) * (damp1_U1[i] + dampt_Ut[i]*w1);
        Rxt_U[i+4] = -0.5 * (xt-x1) * (dampt_Ut[i]*w2);
    }
    Rxt_U[2] -= 1.0;

    Real Ut_x1[4]={0}, Ut_x2[4]={0};
    for (int i = 0; i < 4; ++i) {
        Ut_x1[i] = (U2[i]-U1[i]) * (w2 - 1) / dx;
        Ut_x2[i] = (U2[i]-U1[i]) * (-w2) / dx;
    }
    
    
    // for natural transition the amp is fixed at Ut
    Ut_x1[2] = 0;
    Ut_x2[2] = 0;
    
    
    Real Rxt_x1 = 0.5 * (damp1 + dampt);    
    Real Rxt_x2 = 0.0;                    
    for (int i = 0; i < 4; ++i) {
        Rxt_x1 -= 0.5 * (xt - x1) * dampt_Ut[i] * Ut_x1[i];
        Rxt_x2 -= 0.5 * (xt - x1) * dampt_Ut[i] * Ut_x2[i];
    }

    // sensitivity of xt w.r.t. U,x from Rxt(xt,U,x) = 0 constraint
    Real xt_U[8]={0}, xt_U1[4]={0}, xt_U2[4]={0};
    Real xt_x1 = -Rxt_x1 / Rxt_xt;
    Real xt_x2 = -Rxt_x2 / Rxt_xt;
    
    for (int i = 0; i < 8; ++i){xt_U[i] = -Rxt_U[i] / Rxt_xt;}

    for (int i = 0; i < 4; ++i) {
        xt_U1[i] = xt_U[i];
        xt_U2[i] = xt_U[i+4];
    }

    for (int i = 0; i < 4; ++i) {
        Ut_x1[i] += Ut_xt[i] * xt_x1;
        Ut_x2[i] += Ut_xt[i] * xt_x2;
    }

    // include derivatives w.r.t. xt in Ut_x1 and Ut_x2
    Real Ut_U1[16]={0}, Ut_U2[16]={0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            Ut_U1[i + 4*j] = (i == j ? w1 : 0.0) + (U2[i]-U1[i]) * xt_U1[j] / dx;  // the xt_U1 term goes to zero for fixed
            Ut_U2[i + 4*j] = (i == j ? w2 : 0.0) + (U2[i]-U1[i]) * xt_U2[j] / dx;
        }

    // Copy to Utl, Utt (laminar and turbulent transition states)
    Real Utl[4] = {0}, Utt[4] = {0}, Utl_x1[4] = {0}, Utl_x2[4] = {0}, Utt_x1[4] = {0}, Utt_x2[4] = {0};
    Real Utl_U1[16] = {0}, Utl_U2[16] = {0}, Utt_U1[16] = {0}, Utt_U2[16] = {0};

    
    for (int i = 0; i < 4; ++i) {
        Utl[i] = Ut[i]; Utt[i] = Ut[i];
        Utl_x1[i] = Ut_x1[i]; Utl_x2[i] = Ut_x2[i];
        Utt_x1[i] = Ut_x1[i]; Utt_x2[i] = Ut_x2[i];
        for (int j = 0; j < 4; ++j) {
            Utl_U1[i + 4*j] = Ut_U1[i + 4*j];
            Utl_U2[i + 4*j] = Ut_U2[i + 4*j];
            Utt_U1[i + 4*j] = Ut_U1[i + 4*j];
            Utt_U2[i + 4*j] = Ut_U2[i + 4*j];
        }
    }
    Utl[2] = ncrit;
    Utl_x1[2] = 0; Utl_x2[2] = 0;
    for (int j = 0; j < 4; ++j) {
        Utl_U1[2 + 4*j] = 0.0;
        Utl_U2[2 + 4*j] = 0.0;
    }

    Real cttr, cttr_Ut[4]={0};
    cttr = get_cttr(Ut[0],Ut[1],Ut[2],Ut[3],true,param,cttr_Ut);
    Utt[2] = cttr;
    for (int j = 0; j < 4; ++j) {
        Real sumU1 = 0.0, sumU2 = 0.0;
        for (int i = 0; i < 4; ++i) {
            sumU1 += cttr_Ut[i] * Ut_U1[i + 4*j];
            sumU2 += cttr_Ut[i] * Ut_U2[i + 4*j];
        }
        Utt_U1[2 + 4*j] = sumU1;
        Utt_U2[2 + 4*j] = sumU2;
    }
    Utt_x1[2] = 0.0;
    Utt_x2[2] = 0.0;
    for (int i = 0; i < 4; ++i) {
        Utt_x1[2] += cttr_Ut[i] * Ut_x1[i];
        Utt_x2[2] += cttr_Ut[i] * Ut_x2[i];
    }

    //param.turb = false;
    Real Rl[3]={0}, Rl_U[24]={0}, Rl_x[6]={0};
    residual_station(U1, Utl, x1, xt, aux1, aux2, false, false, false, param, Rl, Rl_U, Rl_x);
    
    //param.turb = true;
    Real Rt[3]={0}, Rt_U[24]={0}, Rt_x[6]={0};
    residual_station(Utt, U2, xt, x2, aux1, aux2, false, true, false, param, Rt, Rt_U, Rt_x);

    for (int i = 0; i < 3; ++i){
        R[i] = Rl[i] + Rt[i];
    }
    // R_U = dRl/dU1 + dRl/dU2 + dRt/dU1 + dRt/dU2
    // Compute contributions using chain rule with intermediate derivatives

    Real Rl_U1[12]={0}, Rl_Utl[12]={0};
    Real Rt_Utt[12]={0}, Rt_U2[12]={0};
    
    // do with pointers ??
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 4; ++c){
            Rl_U1[colMajorIndex(r,c,3)]    = Rl_U[colMajorIndex(r, c, 3)];
            Rl_Utl[colMajorIndex(r,c,3)]   = Rl_U[colMajorIndex(r, c+4, 3)];
            Rt_Utt[colMajorIndex(r,c,3)]   = Rt_U[colMajorIndex(r, c, 3)];
            Rt_U2[colMajorIndex(r,c,3)]    = Rt_U[colMajorIndex(r,c+4,3)];
        }
    }

    // R_x = dRl/dx1 + dRl/dx2 + dRt/dx1 + dRt/dx2
    for (int i = 0; i < 3; ++i) {
        R_x[i] = Rl_x[i] + Rt_x[i] * xt_x1;        // dR/dx1
        R_x[i+3] = Rl_x[i+3] + Rt_x[i] * xt_x2;    // dR/dx2
    }


    Real R_U1[12] = {0},R_U2[12] = {0},tmp1[12]={0};

    // Calculate R_U1
    cnp::add_inplace<12>(R_U1,Rl_U1);
    cnp::matmat_mul<3,4,4>(Rl_Utl, Utl_U1, tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    const Real* Rl_xCol1 = Rl_x + 3;  // Rl_x[:,1]
    cnp::outer_product<3,4>(Rl_xCol1,xt_U1,tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    cnp::matmat_mul<3,4,4>(Rt_Utt,Utt_U1, tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    cnp::outer_product<3,4>(Rt_x,xt_U1,tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    
    // Calculate R_U2
    cnp::add_inplace<12>(R_U2,Rt_U2);
    cnp::matmat_mul<3,4,4>(Rl_Utl,Utl_U2, tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    const Real* Rt_xCol1 = Rt_x + 3;  // Rl_x[:,1]
    cnp::outer_product<3,4>(Rl_xCol1,xt_U2,tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    cnp::matmat_mul<3,4,4>(Rt_Utt,Utt_U2, tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    cnp::outer_product<3,4>(Rt_x,xt_U2,tmp1); cnp::add_inplace<12>(R_U2,tmp1);

    cnp::hstack<12,12>(R_U1,R_U2,R_U); // R_U

    // do R_x: 

    Real R_x1[3]={0},R_x2[3]={0};
    Real tmp2[3]={0};

    cnp::add_inplace<3>(R_x1,Rl_x);
    cnp::scalar_mul<3>(Rl_xCol1,xt_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);
    cnp::scalar_mul<3>(Rt_x,xt_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);
    cnp::matmat_mul<3,4,1>(Rl_Utl,Utl_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);
    cnp::matmat_mul<3,4,1>(Rt_Utt,Utt_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);

    cnp::add_inplace<3>(R_x2,Rt_xCol1);
    cnp::scalar_mul<3>(Rl_xCol1,xt_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);
    cnp::scalar_mul<3>(Rt_x,xt_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);
    cnp::matmat_mul<3,4,1>(Rl_Utl,Utl_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);
    cnp::matmat_mul<3,4,1>(Rt_Utt,Utt_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);

    cnp::hstack<3,3>(R_x1,R_x2,R_x) ;
}


void residual_transition_forced(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Param& param,
    const Real transPos,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6]
) {

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
    
    // gradient of how transition states change wrt starting and ending locations (linear variation across panel)
    Real Ut_x1[4]={0}, Ut_x2[4]={0};
    for (int i = 0; i < 4; ++i) {
        Ut_x1[i] = (U2[i]-U1[i]) * (w2 - 1) / dx;
        Ut_x2[i] = (U2[i]-U1[i]) * (-w2) / dx;
    }
    

    // gradient of transition state wrt start and end states (again linear interp across)
    Real Ut_U1[16]={0}, Ut_U2[16]={0};
    for (int i = 0; i < 4; ++i) {
        int idx = colMajorIndex(i,i,4);  // Diagonal element in column-major
        Ut_U1[idx] = w1;
        Ut_U2[idx] = w2;
    }

    // Copy to Utl, Utt (laminar and turbulent transition states)
    Real Utl[4] = {0}, Utt[4] = {0}, Utl_x1[4] = {0}, Utl_x2[4] = {0}, Utt_x1[4] = {0}, Utt_x2[4] = {0};
    Real Utl_U1[16] = {0}, Utl_U2[16] = {0}, Utt_U1[16] = {0}, Utt_U2[16] = {0};

    // coptying lam/turb states at the transition point
    for (int i = 0; i < 4; ++i) {
        Utl[i] = Ut[i];
        Utt[i] = Ut[i];
        Utl_x1[i] = Ut_x1[i];
        Utl_x2[i] = Ut_x2[i];
        Utt_x1[i] = Ut_x1[i];
        Utt_x2[i] = Ut_x2[i];
        for (int j = 0; j < 4; ++j) {
            Utl_U1[i + 4*j] = Ut_U1[i + 4*j];
            Utl_U2[i + 4*j] = Ut_U2[i + 4*j];
            Utt_U1[i + 4*j] = Ut_U1[i + 4*j];
            Utt_U2[i + 4*j] = Ut_U2[i + 4*j];
        }
    }
    

    // We dont care about the laminar amplicifcation growth at transition as we force it, so remove 
    // residuals for that in Utl
    Utl_x1[2] = 0.0;
    Utl_x2[2] = 0.0;
    for (int j = 0; j < 4; ++j){
        Utl_U1[colMajorIndex(2,j,4)] = 0.0;
        Utl_U2[colMajorIndex(2,j,4)] = 0.0;
    }

    // initialise the turb state with value of shear stress coeff (note Ut[2] has no effect)
    Real cttr, cttr_Ut[4]={0};
    cttr = get_cttr(Ut[0],Ut[1],Ut[2],Ut[3],true,param,cttr_Ut);
    Utt[2] = cttr;
    

    // doing that dUtt/dU1 [2,:] = dot(dcttr/dUt, dUt/dU1)
    // doing that dUtt/dU2 [2,:] = dot(dcttr/dUt, dUt/dU2)
    // but dUt/dU1 is just diaginal with value of w1 so code is simpler: 
    for (int j = 0; j < 4; ++j) {
        Utt_U1[colMajorIndex(2,j,4)] = cttr_Ut[j] * w1;
        Utt_U2[colMajorIndex(2,j,4)] = cttr_Ut[j] * w2;
    }
    
    

    // using dUtt/dx1 [2,:] = dot(dcttr/dUt, dUt/dx1)
    // using dUtt/dx2 [2,:] = dot(dcttr/dUt, dUt/dx2)
    Utt_x1[2] = 0.0;
    Utt_x2[2] = 0.0;
    for (int i = 0; i < 4; ++i) {
        Utt_x1[2] += cttr_Ut[i] * Ut_x1[i];
        Utt_x2[2] += cttr_Ut[i] * Ut_x2[i];
    }

    // residual for laminar section (remember that dont care about amplifcation residual here)
    Real Rl[3]={0}, Rl_U[24]={0}, Rl_x[6]={0};
    residual_station(U1, Utl, x1, xt, 0.0, 0.0, false, false, false, param, Rl, Rl_U, Rl_x);
    
    Rl[2] = 0.0;
    for (int i = 0; i < 8; ++i) {
        Rl_U[colMajorIndex(2,i,3)] = 0.0;
    }
     for (int i = 0; i < 2; ++i) {
        Rl_x[colMajorIndex(2,i,3)] = 0.0;
    }

    // residual for the turbulent section, here residual for shear stress matters
    Real Rt[3]={0}, Rt_U[24]={0}, Rt_x[6]={0};
    residual_station(Utt, U2, xt, x2, 0.0, 0.0, false, true, false, param, Rt, Rt_U, Rt_x);

    // sum residuals for total panel residuals
    for (int i = 0; i < 3; ++i){
        R[i] = Rl[i] + Rt[i];
    }

    Real Rl_U1[12]={0}, Rl_Utl[12]={0};
    Real Rt_Utt[12]={0}, Rt_U2[12]={0};
    
    // seperating Rl_U correctly and Rt_U
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 4; ++c){
            Rl_U1[colMajorIndex(r,c,3)]    = Rl_U[colMajorIndex(r, c, 3)];
            Rl_Utl[colMajorIndex(r,c,3)]   = Rl_U[colMajorIndex(r, c+4, 3)];
            Rt_Utt[colMajorIndex(r,c,3)]   = Rt_U[colMajorIndex(r, c, 3)];
            Rt_U2[colMajorIndex(r,c,3)]    = Rt_U[colMajorIndex(r,c+4,3)];
        }
    }

    // R_x = dRl/dx1 + dRl/dx2 + dRt/dx1 + dRt/dx2
    for (int i = 0; i < 3; ++i) {
        R_x[i] = Rl_x[i];       // dR/dx1
        R_x[i+3] = Rl_x[i+3];  // dR/dx2
    }


    Real R_U1[12] = {0},R_U2[12] = {0},tmp1[12]={0};

    // Calculate R_U1
    cnp::add_inplace<12>(R_U1,Rl_U1);
    cnp::matmat_mul<3,4,4>(Rl_Utl, Utl_U1, tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    cnp::matmat_mul<3,4,4>(Rt_Utt,Utt_U1, tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    
    // Calculate R_U2
    cnp::add_inplace<12>(R_U2,Rt_U2);
    cnp::matmat_mul<3,4,4>(Rl_Utl,Utl_U2, tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    cnp::matmat_mul<3,4,4>(Rt_Utt,Utt_U2, tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    
    cnp::hstack<12,12>(R_U1,R_U2,R_U); // combine into single Residual jacobian

    // do R_x: 
    Real R_x1[3]={0},R_x2[3]={0};
    Real tmp2[3]={0};

    const Real* Rl_xCol1 = Rl_x + 3;  // Rl_x[:,1]
    const Real* Rt_xCol1 = Rt_x + 3;  // Rt_x[:,1]

    cnp::add_inplace<3>(R_x1,Rl_x); // add Rl_x[:,0]
    cnp::matmat_mul<3,4,1>(Rl_Utl,Utl_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);
    cnp::matmat_mul<3,4,1>(Rt_Utt,Utt_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);

    cnp::add_inplace<3>(R_x2,Rt_xCol1);
    cnp::matmat_mul<3,4,1>(Rl_Utl,Utl_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);
    cnp::matmat_mul<3,4,1>(Rt_Utt,Utt_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);

    cnp::hstack<3,3>(R_x1,R_x2,R_x) ;
}