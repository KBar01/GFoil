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
template<typename Real> struct Vsol;

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
template<typename Real> struct Foil;
template<typename Real> struct Param;
template<typename Real> struct Post;
template<typename Real> struct Oper;
template<typename Real> struct Geom;
template<typename Real> struct Wake;
template<typename Real> struct Glob;
template<typename Real> struct Trans;




// init_thermo moved to src/include/solver_funcs.hpp

// space_wake_nodes moved to src/include/solver_funcs.hpp


template<typename Real>
void build_wake(const Foil<Real>& foil, const Geom<Real>& geom, const Oper<Real>& op, Isolc<Real>& isol, Wake<Real>& wake){


    Real firstPanelSize = 0.5*(foil.s[1]-foil.s[0] + foil.s[Ncoords-1]-foil.s[Ncoords-2]);
    Real wakeLength = geom.wakelen*geom.chord ;
    Real wakePanelSizes[Nwake] ;
    space_wake_nodes(wakeLength,firstPanelSize,wakePanelSizes,foil,wake); //fills wake panel sizes

    // dont need to make xyw or twm use that in wake Struct through modification
    Real midpointTE[2];
    midpointTE[0] = 0.5*(foil.x[colMajorIndex(0,0,2)] + foil.x[colMajorIndex(0,Ncoords-1,2)]);
    midpointTE[1] = 0.5*(foil.x[colMajorIndex(1,0,2)] + foil.x[colMajorIndex(1,Ncoords-1,2)]);
    
    Real normTE[2];
    Real tangTE[2];
    normTE[0] = foil.x[colMajorIndex(0,Ncoords-1,2)] - foil.x[colMajorIndex(0,0,2)];
    normTE[1] = foil.x[colMajorIndex(1,Ncoords-1,2)] - foil.x[colMajorIndex(1,0,2)];

    tangTE[0] = normTE[1];
    tangTE[1] = -normTE[0];

    // Fill first wake coordinates
    wake.x[0] = midpointTE[0] + (1.0e-5)*tangTE[0]*geom.chord ;
    wake.x[1] = midpointTE[1] + (1.0e-5)*tangTE[1]*geom.chord ;

    
    
    Real v1[2],v2[2],v1Store[2];
    Real norm,wakeLengthsDiff ; 
    for (int i=0; i<Nwake-1; ++i){
        
        v1Store[0]=0.0;
        v1Store[1]=0.0;
        v2[0]=0.0;
        v2[1]=0.0;
        // Consider storing the results of these for uewi
        inviscid_velocity(foil,isol.gammas,op.Vinf,op.alpha,wake.x[colMajorIndex(0,i,2)],wake.x[colMajorIndex(1,i,2)],v1Store);
        
        norm = norm2(v1Store);
        v1[0] = v1Store[0]/norm;
        v1[1] = v1Store[1]/norm;
        wake.t[colMajorIndex(0,i,2)] = v1[0];
        wake.t[colMajorIndex(1,i,2)] = v1[1];
        
        // Store inviscid edge velocity of wake point
        isol.uewi[i] = (v1Store[0] * v1[0]) + (v1Store[1] * v1[1]);

        //Predictor step:
        wakeLengthsDiff = wakePanelSizes[i+1]-wakePanelSizes[i];
        wake.x[colMajorIndex(0,i+1,2)] = wake.x[colMajorIndex(0,i,2)] + wakeLengthsDiff*v1[0];
        wake.x[colMajorIndex(1,i+1,2)] = wake.x[colMajorIndex(1,i,2)] + wakeLengthsDiff*v1[1];

        inviscid_velocity(foil,isol.gammas,op.Vinf,op.alpha,wake.x[colMajorIndex(0,i+1,2)],wake.x[colMajorIndex(1,i+1,2)],v2);
        norm = norm2(v2);
        v2[0] /= norm;
        v2[1] /= norm;
        wake.t[colMajorIndex(0,i+1,2)] = v2[0];
        wake.t[colMajorIndex(1,i+1,2)] = v2[1];
        
        //Corrector step:
        wake.x[colMajorIndex(0,i+1,2)] = wake.x[colMajorIndex(0,i,2)] + wakeLengthsDiff*0.5*(v1[0]+v2[0]);
        wake.x[colMajorIndex(1,i+1,2)] = wake.x[colMajorIndex(1,i,2)] + wakeLengthsDiff*0.5*(v1[1]+v2[1]);
    }

    // Fill last wake node not done in loop
    v1Store[0]=0.0;
    v1Store[1]=0.0;
    inviscid_velocity(foil,isol.gammas,op.Vinf,op.alpha,wake.x[colMajorIndex(0,Nwake-1,2)],wake.x[colMajorIndex(1,Nwake-1,2)],v1Store);
    norm = norm2(v1Store);
    v1[0] = v1Store[0]/norm;
    v1[1] = v1Store[1]/norm;
    wake.t[colMajorIndex(0,Nwake-1,2)] = v1[0];
    wake.t[colMajorIndex(1,Nwake-1,2)] = v1[1];
    
    isol.uewi[Nwake-1] = (v1Store[0] * v1[0]) + (v1Store[1] * v1[1]);

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
void finishdRdU_AD(const Foil<Real>&foil, const Isolc<Real>&isolc, const Isolv<Real>&isolv, Glob<Real>& glob, Vsol<Real>& vsol, const Oper<Real>& oper) {
    
    
    constexpr int Nsys = Ncoords+Nwake;

    // Step 1: Modify ue array to avoid 0 or negative
    int nrows = 4; // Since U is shaped (4, Nsys) in column-major
    Real ue[Nsys] = {0};
    Real uemax = 0.0;
    for (int i = 0; i < Nsys; ++i){
        uemax = std::max(uemax, std::abs(glob.U[colMajorIndex(3,i,4)]));
    }
    for (int i = 0; i < Nsys; ++i){
        ue[i] = std::max(glob.U[colMajorIndex(3,i,4)], 1e-10*uemax);
    }

    // Step 2: Get ueinv
    Real ueinv[Nsys]={0};
    get_ueinv(isolc,isolv,ueinv);

    // Step 3: Build Residual R
    Real ds[Nsys];
    for (int i = 0; i < Nsys; ++i) {ds[i] = glob.U[colMajorIndex(1, i, 4)];}

    Real tempRHS[Nsys];
    cnp::mul<Nsys>(ds,ue,tempRHS); // ds*ue

    Real* Rpointer = &glob.R[3*Nsys] ; 
    cnp::matmat_mul<Nsys,Nsys,1>(vsol.ue_m,tempRHS,Rpointer);

    for (int i = 0; i < Nsys; ++i){
        Rpointer[i] = ue[i] - (ueinv[i] + Rpointer[i]);
    }
}

template<typename Real>
void stagpoint_move_AD(Isolv<Real>& isol,Glob<Real>& glob,const Foil<Real>& foil,const Wake<Real>& wake,Vsol<Real>&vsol,const int (&currStag)[2]) {
    
    const int* I = currStag;       // pointer to stagnation indices (keeps it cleaner)
    bool newpanel = true;
    
    // Compute new stagnation location
    Real u0 = glob.U[colMajorIndex(3,I[0],4)]; // velocities at stag point panel
    Real u1 = glob.U[colMajorIndex(3,I[1],4)]; 

    Real den = u0 + u1;
    Real w1 = u1 / den;
    Real w2 = u0 / den;

    isol.stagArcLocation = w1*foil.s[I[0]] + w2*foil.s[I[1]]; // new arclength location

    for (int d=0; d<2; ++d) {   // new x coord of stagnation point
        isol.stagXLocation[d] = w1*foil.x[colMajorIndex(d,I[0],2)] + w2*foil.x[colMajorIndex(d,I[1],2)];
    }

    Real ds = foil.s[I[1]] - foil.s[I[0]];
    isol.sstag_ue[0] = u1 * ds / (den * den);
    isol.sstag_ue[1] = -u0 * ds / (den * den);


    // updating the array of arclength from stagnation at every node


    cnp::scalar_sub_abs<Ncoords>(foil.s,isol.stagArcLocation,isol.distFromStag);
    Real* xiWake = isol.distFromStag + Ncoords ;
    cnp::scalar_sub<Nwake>(wake.s,isol.stagArcLocation,xiWake);

    for (int i=0; i<=I[0]; ++i){ isol.edgeVelSign[i] = -1;}
    for (int i=I[0]+1; i<Ncoords; ++i){ isol.edgeVelSign[i] = 1;}

    isol.stagIndex[0] = I[0] ;
    isol.stagIndex[1] = I[1] ;


    identify_surfaces(isol,vsol);
    rebuild_ue_m(foil,wake,isol,vsol,false);
    
};

template<typename Real>
void stagnation_state(const Real*U1,const Real*U2,const Real x1,const Real x2,
    Real (&Ust)[4],Real&xst){


    Real dx = x2-x1;
    Real dx_x[2] = {-1, 1};
    Real rx = x2/x1;
    Real rx_x[2] = {-rx/x1,1/x1};
  
    // linear extrapolation weights and stagnation state
    Real w1 =  x2/dx, w1_x[2] = {-w1/dx*dx_x[0], -w1/dx*dx_x[1] + 1/dx};
    Real w2 = -x1/dx, w2_x[2] = {-w2/dx*dx_x[0] -1/dx, -w2/dx*dx_x[1]};
        
    for (int i=0;i<4;++i){Ust[i] = U1[i]*w1 + U2[i]*w2;}
    // quadratic extrapolation of the edge velocity for better slope, ue=K*x
    Real wk1 = rx/dx,       wk1_x[2] = {rx_x[0]/dx - wk1/dx*dx_x[0], rx_x[1]/dx - wk1/dx*dx_x[1]};
    Real wk2 = -1/(rx*dx),  wk2_x[2] = {-wk2*(rx_x[0]/rx + dx_x[0]/dx), -wk2*(rx_x[1]/rx + dx_x[1]/dx)}; 
    Real K = wk1*U1[3] + wk2*U2[3] ;
  
    //stagnation coord cannot be zero, but must be small
    xst = 1e-6;
    Ust[3] = K*xst ; // linear dep of ue on x near stagnation

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