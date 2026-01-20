#pragma once

#include "real_type.hpp"
#include "data_structs.hpp"
#include "vector_ops.hpp"
#include "get_funcs.hpp"
#include "panel_funcs.hpp"
#include "residuals.hpp"


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




template<typename Real>
void init_thermo(const Oper<Real>& oper,Param<Real>& param,const Geom<Real>& geom){

    param.Vinf = oper.Vinf;
    param.muinf = oper.rho*oper.Vinf*geom.chord/oper.Re;
    if (oper.Ma > 0.0) {
        // Compressibility corrections
        Real gmi = param.gam - 1.0;
        param.KTb = std::sqrt(1.0 - oper.Ma*oper.Ma);
        param.KTl = (oper.Ma * oper.Ma) / std::pow(1.0 + param.KTb, 2);
        param.H0 = ((1.0 + 0.5 * gmi*oper.Ma*oper.Ma)*oper.Vinf*oper.Vinf) / (gmi*oper.Ma*oper.Ma);
        
        Real Tr = 1.0 - (0.5 * oper.Vinf*oper.Vinf) / param.H0;  // Freestream/Stagnation temperature ratio
        Real finf = std::pow(Tr, 1.5) * (1.0 + param.Tsrat) / (Tr + param.Tsrat);  // Sutherland's ratio
        
        param.cps = (2.0 / (param.gam * oper.Ma*oper.Ma)) * 
                    (std::pow((1.0 + 0.5*gmi*oper.Ma*oper.Ma) / (1.0 + 0.5 * gmi), param.gam / gmi) - 1.0);
        
        param.mu0 = param.muinf / finf;  // Stagnation viscosity
        param.rho0 = oper.rho * std::pow(1.0 + 0.5 * gmi * oper.Ma * oper.Ma, 1.0 / gmi);  // Stagnation density
    } else {
        // Incompressible case
        param.mu0 = param.muinf;
        param.rho0 = oper.rho;
    }
}

template<typename Real>
void space_wake_nodes(const Real& wakeLength,const Real& firstPanelLength,Real* wakeSpacing,const Foil<Real>&foil,Wake<Real>&wake){

    int Nintervals = Nwake -1;
    Real d = wakeLength/firstPanelLength ;
    Real a = Nintervals*(Nintervals-1.0)*(Nintervals-2.0)/6.0 ;
    Real b = Nintervals*(Nintervals-1.0)/2.0 ;
    Real c = Nintervals-d;

    Real disc = std::max(b*b-4.0*a*c, 0.0);
    Real r    = 1 + (-b + std::sqrt(disc))/(2*a) ;

    // Newton-Raphson iterations
    Real R, R_r, dr;
    for (int k = 0; k < 10; ++k) {
        
        R = std::pow(r, Nintervals) - 1 - d * (r - 1);
        R_r = Nintervals * std::pow(r, Nintervals - 1) - d;
        dr = -R / R_r;
        if (std::abs(dr) < 1e-6) break;
        r -= R / R_r;
    }


    wakeSpacing[0] = 0.0 ;
    Real foilEndS = foil.s[Ncoords-1];
    wake.s[0] = foilEndS + 0.0;

    // Compute cumulative sum
    Real term = firstPanelLength;

    for (int i = 0; i < Nwake-1; ++i) {
        wakeSpacing[i+1] = wakeSpacing[i] + term;
        wake.s[i+1] = foilEndS + wakeSpacing[i+1];
        // include r multiplication to avoid constantly calling power
        term *=r;
    }
}


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
void stagpoint_find(const Isolc<Real>& isolc, Isolv<Real> & isolv, const Foil<Real>& foil,const Wake<Real>& wake) {


    int j=0;
    // Find first positive gamma
    for (; j < Ncoords; ++j) {
        if (isolc.gammas[j] > 0) break;
    }
    int I[2] = {j-1, j};
    isolv.stagIndex[0] = I[0];
    isolv.stagIndex[1] = I[1];

    Real G[2] = {isolc.gammas[I[0]], isolc.gammas[I[1]]};
    Real S[2] = {foil.s[I[0]], foil.s[I[1]]};

    Real den = G[1] - G[0];
    Real w1 = G[1] / den;
    Real w2 = -G[0] / den;

    isolv.stagArcLocation = w1 * S[0] + w2 * S[1];

    // Compute x location in column-major form
    isolv.stagXLocation[0] = foil.x[colMajorIndex(0,j-1,2)] * w1 + foil.x[colMajorIndex(1,j,2)] * w2;   // x-coordinate
    isolv.stagXLocation[1] = foil.x[colMajorIndex(1,j-1,2)] * w1 + foil.x[colMajorIndex(1,j,2)] * w2; // y-coordinate

    // Compute sign conversion (sgnue)
    for (int i = j; i < Ncoords; ++i) {
        isolv.edgeVelSign[i] = 1.0;
    }

    // Compute distance from stagnation point (xi)
    for (int i = 0; i < Ncoords; ++i) {
        isolv.distFromStag[i] = std::abs(foil.s[i] - isolv.stagArcLocation);
    }

    for (int i = 0; i < Nwake; ++i) {
        isolv.distFromStag[Ncoords + i] = wake.s[i] - isolv.stagArcLocation;
    }

}

// Helper function to mimic Python's range(start, end, step)
std::vector<int> range(int start, int end, int step = 1) {
    std::vector<int> result;
    if (step > 0) {
        for (int i = start; i < end; i += step)
            result.push_back(i);
    } else if (step < 0) {
        for (int i = start; i > end; i += step)
            result.push_back(i);
    }
    return result;
}


// Core function to fill in surface indices
template<typename Real>
void identify_surfaces(const Isolv<Real>& isolv, Vsol<Real>& vsol) {
    vsol.Is.clear(); // Clear any previous data

    // Lower surface (reverse order from Istag[0] to 0)
    vsol.Is.push_back(range(isolv.stagIndex[0], -1, -1));

    // Upper surface (from Istag[1] to foil.N-1)
    vsol.Is.push_back(range(isolv.stagIndex[1], Ncoords));

    // Wake (from foil.N to foil.N + wake.N - 1)
    vsol.Is.push_back(range(Ncoords, Ncoords+Nwake));
}

template<typename Real>
void set_wake_gap(const Foil<Real>&foil,const Isolv<Real>&isol,Vsol<Real>&vsol){

    Real lengthScaleFactr = 2.5;
    Real dtdx = std::min(std::max(foil.te.dtdx,-3.0/lengthScaleFactr),3.0/lengthScaleFactr);
    Real Lw =lengthScaleFactr*foil.te.hTE;

    Real xib ;
    for (int i=0; i<Nwake;++i){
        xib = (isol.distFromStag[Ncoords+i] - isol.distFromStag[Ncoords])/Lw;
        if (xib <= 1.0){
            vsol.wgap[i] = foil.te.hTE*(1+(2+lengthScaleFactr*dtdx)*xib)*(1-xib)*(1-xib) ;
        }
    }
}


template<typename Real>
void calc_force(const Oper<Real>&op, const Geom<Real>&geom, const Param<Real>&par,const Foil<Real>&foil, const Glob<Real>&glob, Post<Real>& post) {

    // Compute dynamic pressure
    Real qinf = 0.5*op.rho *(op.Vinf * op.Vinf);
    constexpr int Nsys = Ncoords+Nwake;
    
    // Calculate the pressure coefficient at each node
    Real ue[Nsys]={0};
    for (int i=0;i<Nsys;++i){ue[i] = glob.U[colMajorIndex(3,i,4)];}
    get_cp(post,op,par,ue);     // assigns cp in post struct

    // Pre-allocate variables for memory efficiency
    Real cos_alpha = std::cos(op.alpha);
    Real sin_alpha = std::sin(op.alpha);
    
    // Loop over panels for integration
    for (int i0 = 1; i0 <= Ncoords; ++i0) {
        
        int i = (i0 == Ncoords) ? 0 : i0;
        int ip = (i0 == Ncoords) ? Ncoords-1 : i0-1;

        const Real x1[2] = {foil.x[colMajorIndex(0,ip,2)],foil.x[colMajorIndex(1,ip,2)]}; // points
        const Real x2[2] = {foil.x[colMajorIndex(0,i ,2)],foil.x[colMajorIndex(1,i ,2)]};

        Real dxv[2] = {x2[0] - x1[0], x2[1] - x1[1]};
        Real dx1[2] = {x1[0] - geom.xref[0], x1[1] - geom.xref[1]};
        Real dx2[2] = {x2[0] - geom.xref[0], x2[1] - geom.xref[1]};

        Real dx1nds = dxv[0] * dx1[0] + dxv[1] * dx1[1];
        Real dx2nds = dxv[0] * dx2[0] + dxv[1] * dx2[1];
        //Real dx1nds = (x1[0] - geom.xref[0]) * dxv[1] + (x1[1] - geom.xref[1]) * dxv[0];
        //Real dx2nds = (x2[0] - geom.xref[0]) * dxv[1] - (x2[1] - geom.xref[1]) * dxv[0];

        Real dx = -dxv[0] * cos_alpha - dxv[1] * sin_alpha;
        Real dz =  dxv[1] * cos_alpha - dxv[0] * sin_alpha;

        Real cp1 = post.cp[ip], cp2 = post.cp[i];
        Real cpbar = 0.5 * (cp1 + cp2);

        post.cl += dx * cpbar;
        //post.cl_ue[ip] += dx * 0.5 * cp_ue[ip];
        //post.cl_ue[i]  += dx * 0.5 * cp_ue[i];

        //cl_alpha += cpbar * (sind(alpha) * dxv[0] - cosd(alpha) * dxv[1]) * deg2rad;

        post.cm += cp1 * dx1nds / 3.0 + cp1 * dx2nds / 6.0 + cp2 * dx1nds / 6.0 + cp2 * dx2nds / 3.0;
        //cdpi += dz * cpbar;
    }

    // Normalize by chord
    post.cl /= geom.chord;
    post.cm /= (geom.chord*geom.chord);
    
    int iw = Ncoords+Nwake-1;  // end of wake
    const Real* U = &glob.U[colMajorIndex(0,iw,4)];

    Real H, H_U[4];
    H = get_H(U[0],U[1],H_U);

    Real uk, uk_ue;
    uk = get_uk(U[3],par,uk_ue);

    post.cd = 2.0 * U[0] * pow(uk / op.Vinf, (5.0 + H) / 2.0);
    
}

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
    cnp::mul<Real,Nsys>(ds,ue,tempRHS); // ds*ue

    Real* Rpointer = &glob.R[3*Nsys] ; 
    cnp::matmat_mul<Real,Nsys,Nsys,1>(vsol.ue_m,tempRHS,Rpointer);

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


    cnp::scalar_sub_abs<Real,Ncoords>(foil.s,isol.stagArcLocation,isol.distFromStag);
    Real* xiWake = isol.distFromStag + Ncoords ;
    cnp::scalar_sub<Real,Nwake>(wake.s,isol.stagArcLocation,xiWake);

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


