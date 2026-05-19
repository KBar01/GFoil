#pragma once

#include "real_type.hpp"
#include "data_structs.hpp"
#include "vector_ops.hpp"
#include "get_funcs.hpp"
#include "panel_funcs.hpp"
#include "residuals.hpp"
#include "calc_ue_m.hpp"  // shared calc_ue_m template; also provides solve_sys_ue
#include "solver_funcs.hpp"  // shared duck-typed template implementations


template<typename Real> struct Isolc;
template<typename Real> struct Isolv;
// Vsol (and other types) are now template aliases — no struct forward-decl.

// AD version of get_ueinv: kept here because it takes Isolc<Real>+Isolv<Real>
// (the Isol/Isolc split means this cannot be in the shared src/include/get_funcs.hpp).
template<typename Real>
void get_ueinv(const Isolc<Real>& isolc, const Isolv<Real>& isolv, Real* ueinv) {
    for (int i = 0; i < Ncoords; ++i) {
        ueinv[i] = isolv.edgeVelSign[i] * isolc.gammas[i];
    }
    for (int i = 0; i < Nwake; ++i) {
        ueinv[Ncoords+i] = isolc.uewi[i];
    }
    ueinv[Ncoords] = ueinv[Ncoords-1];
}
// These are now template aliases in data_structs.hpp — no struct forward-decls needed.




// init_thermo moved to src/include/solver_funcs.hpp

// space_wake_nodes moved to src/include/solver_funcs.hpp


template<typename Real>
void build_wake(const Foil<Real>& foil, const Geom<Real>& geom,
                const Oper<Real>& op, Isolc<Real>& isol, Wake<Real>& wake) {
    build_wake_impl<>(foil, geom, op, isol, wake);
}

template<typename Real>
void stagpoint_find(const Isolc<Real>& isolc, Isolv<Real>& isolv,
                    const Foil<Real>& foil, const Wake<Real>& wake) {
    // isolc = gamma source, isolv = variable destination; no sstag_g in AD
    stagpoint_find_impl<false>(isolc, isolv, foil, wake);
}

// range() and identify_surfaces moved to src/include/solver_funcs.hpp

// set_wake_gap moved to src/include/solver_funcs.hpp


// calc_force moved to src/include/solver_funcs.hpp

template<typename Real>
void finishdRdU_AD(const Foil<Real>&foil, const Isolc<Real>&isolc,
                   const Isolv<Real>&isolv, Glob<Real>& glob,
                   Vsol<Real>& vsol, const Oper<Real>& oper) {
    ue_residual_kernel<Real>(isolc, isolv, vsol, glob);
}

template<typename Real>
void stagpoint_move_AD(Isolv<Real>& isol, Glob<Real>& glob,
                       const Foil<Real>& foil, const Wake<Real>& wake,
                       Vsol<Real>& vsol, const int (&currStag)[2]) {
    stagpoint_move_impl<Real>(isol, glob, foil, wake, vsol, currStag);
}

template<typename Real>
void stagnation_state(const Real*U1, const Real*U2, const Real x1, const Real x2,
    Real (&Ust)[4], Real& xst) {
    stagnation_state_impl<Real>(U1, U2, x1, x2, Ust, xst);
}

template<typename Real>
void build_glob_RV_AD(const Foil<Real>&foil, const Vsol<Real>&vsol,const Isolv<Real>&isol,Glob<Real>&glob, Param<Real>&param, Trans<Real>&tdata){
    
    constexpr int RVsize = 4*(Ncoords+Nwake);
    constexpr int RXsize = 3*(Ncoords+Nwake);
    
    const Real* xi = isol.distFromStag;
    for (int si = 0; si < 3; ++si) {    // for each surface (upper/lower/wake)
        

        const std::vector<int>& Is = vsol.Is[si]; // list of surface node indices from stag point
        const int nSurfPoints = Is.size();

        // Check for edge case of first node hitting stag point exactly
        // i0 will be 1 if this happens
        int i0 = ((si < 2) && (xi[Is[0]] < 1e-8 * xi[Is[nSurfPoints-1]])) ? 1 : 0;

        bool turb = false, wake=false ; 
        
        if (si < 2) {
            
            Real R1[3];
            
            Real* U1 = &glob.U[colMajorIndex(0,Is[i0],4)];
            Real* U2 = &glob.U[colMajorIndex(0,Is[i0+1],4)];
            
            // Compute stagnation state
            Real Ust[4], xst; 
            stagnation_state(U1,U2,xi[Is[i0]],xi[Is[i0+1]],Ust,xst);
            
            Real Ust_copy[4];

            for (int i=0;i<4;++i){
                Ust_copy[i] = Ust[i];
            }
            residual_station<Real>(Ust,Ust_copy,xst,xst,0,0,false,false,true,param,R1);

            
            int J[2] = {Is[i0],Is[i0+1]};

            int Ig = 3*Is[i0];
            for (int row=0;row<3;++row){
                glob.R[Ig+row] = R1[row];
            }
        } 
        else {  // dealing with start of wake
            Real R1[3];
            wake_sys(vsol,foil,glob,param,R1);

            int Ig = 3*Is[i0];
            for (int row=0;row<3;++row){glob.R[Ig+row] = R1[row];}

            wake = true;
            turb = true;
        }

        // loop over remaining points in surface
        for (int i = i0 + 1; i < nSurfPoints; ++i) {
            
            int prevI = i - 1;
            int currI = i;
            bool tran = vsol.turb[Is[prevI]] ^ vsol.turb[Is[currI]];
            
            Real Ri[3];
            
            Real* Uprev = &glob.U[colMajorIndex(0,Is[prevI],4)];
            Real* Ucurr = &glob.U[colMajorIndex(0,Is[currI],4)];

            if (tran){

                int isForced = tdata.isForced[si];
                Real transPos = tdata.transPos[si];

                if (!isForced){
                    residual_transition<Real>(Uprev,Ucurr,xi[Is[prevI]],xi[Is[currI]],0,0,param,Ri);
                }
                else {

                    // do linear interp to find y value
                    Real x1 = foil.x[colMajorIndex(0,Is[prevI],2)],y1 = foil.x[colMajorIndex(1,Is[prevI],2)];
                    Real x2 = foil.x[colMajorIndex(0,Is[currI],2)],y2 = foil.x[colMajorIndex(1,Is[currI],2)];

                    Real yt = y1 + ((tdata.transPos[si] - x1) / (x2 - x1)) * (y2 - y1);                
                    Real distFromStagTrans =  isol.distFromStag[Is[prevI]] + std::sqrt((tdata.transPos[si]-x1)*(tdata.transPos[si]-x1) + (yt-y1)*(yt-y1));

                    residual_transition_forced(Uprev,Ucurr,xi[Is[prevI]],xi[Is[currI]],param,distFromStagTrans,Ri);
                }
            }
            else {
                Real aux1=0,aux2=0;
                if (wake){aux1=vsol.wgap[Is[prevI]-Ncoords];aux2 = vsol.wgap[Is[currI]-Ncoords] ;}
                residual_station(Uprev,Ucurr,xi[Is[prevI]],xi[Is[currI]],aux1,aux2,wake,turb,false,param,Ri);
            }
            
            // update residuals
            int Ig = 3*Is[i];
            for (int j = 0; j < 3; ++j){
                glob.R[Ig + j] += Ri[j];
            }

            if (tran) {turb = true;}
        }
    }
}



// rebuild_ue_m moved to src/include/solver_funcs.hpp