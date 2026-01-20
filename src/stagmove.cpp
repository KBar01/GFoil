#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "vector_ops.hpp"
#include "main_func.h"
#include "data_structs.h"


void stagpoint_move(Isol& isol,Glob& glob,const Foil& foil,const Wake& wake,Vsol&vsol) {
    
    

    int* I = isol.stagIndex;       // pointer to stagnation indices (keeps it cleaner)
    bool newpanel = true;
    if (glob.U[colMajorIndex(3,I[1],4)] < 0) {

        // condition for moving stagnation up a panel, as velocity direction on first top surface panel
        // has reversed (negative sign)
        int j;
        for (j = I[1]; j < Ncoords; ++j) {
            if (glob.U[colMajorIndex(3,j,4)] > 0) break; // find where new top surface begins from ue sign
        }

        int I1 = j; // new start of top surface 
        for (j = I[1]; j < I1; ++j){ glob.U[colMajorIndex(3,j,4)] *= -1.0;} // set everything in changed range back to pos ue
        I[0] = I1 - 1;
        I[1] = I1;          // update the stagnation point in isol struct
    } 
    else if (glob.U[colMajorIndex(3,I[0],4)] < 0) {
        
        // This is moving stagnation down, so further along bottom surface, tell by sign swap again
        int j;
        for (j = I[0]; j >= 0; --j) {
            if (glob.U[colMajorIndex(3,j,4)] > 0) break; // going backwards, find where new bottom surface starts
        }
        assert(j > 0 && "no stagnation point");
        int I0 = j;
        for (j = I0+1; j <= I[0]; ++j) {glob.U[colMajorIndex(3,j,4)] *= -1.0; }
        I[0] = I0;      // set new stagnation in isol struct
        I[1] = I0 + 1;
    } else {
        newpanel = false;   // no change detected, keep current stag point
    }

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


    // If moved to a new panel, update ue signs and recompute surfaces
    if (newpanel) {

        for (int i=0; i<=I[0]; ++i){ isol.edgeVelSign[i] = -1;}
        for (int i=I[0]+1; i<Ncoords; ++i){ isol.edgeVelSign[i] = 1;}

        identify_surfaces(isol,vsol);
        rebuild_ue_m(foil,wake,isol,vsol,true);
    }
}

