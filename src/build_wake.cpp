#include <iostream>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"


void space_wake_nodes(const Real& wakeLength,const Real& firstPanelLength,Real* wakeSpacing,const Foil&foil,Wake&wake){

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



void build_wake(const Foil& foil, const Geom& geom, const Oper& op, Isol& isol, Wake& wake){


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


void set_wake_gap(const Foil&foil,const Isol&isol,Vsol&vsol){

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