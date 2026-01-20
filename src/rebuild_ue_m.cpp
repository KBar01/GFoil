#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"





Real sigma_m[2*(Ncoords-1)];

void rebuild_ue_m(const Foil&foil,const Wake&wake,const Isol&isol,Vsol&vsol,bool realloc){


    // Functions builds ue_m from ue_sigma, sigma_m and edgevelocity signs.
    // This looks more complex as taking advantage of known sparsity pattern
    // of sigma_m, therefore removing unnessessary computation.

    // TODO: mayb make iterators in loops more clear etc
    
    if (realloc){
        for (int i=0;i<(Ncoords+Nwake)*(Ncoords+Nwake);++i){
            vsol.ue_m[i] = 0;
        }
    }

    // Creation of sigma_m, holding non-zero vals only
    for (int i=0; i<Ncoords-1;++i){

        Real ds = foil.s[i+1] - foil.s[i];
        int ind = 2*i;
        sigma_m[ind] = -1.0*isol.edgeVelSign[i] / ds;
        sigma_m[ind+1] = 1.0*isol.edgeVelSign[i+1] / ds;
    }

    // first column in ue_m
    for (int row=0; row<Ncoords+Nwake;++row){
        vsol.ue_m[row] = vsol.ue_sigma[row]*sigma_m[0];
    }

    // columns 1 to Ncoords-1
    int ue_mInt = Ncoords+Nwake;
    int sigma_mIndex = 1;
    for (int col=1; col<=Ncoords-2;++col){

        Real sigma_mValue = sigma_m[sigma_mIndex];

        for (int row=0; row<Ncoords+Nwake; ++row){
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,col-1,Ncoords+Nwake)]*sigma_mValue;
        }

        sigma_mValue = sigma_m[sigma_mIndex+1];
        for (int row=0; row<Ncoords+Nwake; ++row){
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,col,Ncoords+Nwake)]*sigma_mValue;
        }

        sigma_mIndex += 2;
        ue_mInt += (Ncoords+Nwake);
    }

    // Ncoords-1 column in ue_m (0 index)
    int start = (Ncoords+Nwake)*(Ncoords-1);
    for (int row=0; row<Ncoords+Nwake;++row){
        vsol.ue_m[colMajorIndex(row,Ncoords-1,Ncoords+Nwake)] = vsol.ue_sigma[colMajorIndex(row,Ncoords-2,Ncoords+Nwake)]*sigma_m[2*(Ncoords-1)-1];
    }


    // Now contributions from wake section

    Real sigma_mWake[2*(Nwake-1)];

    for (int i=0; i<Nwake-1;++i){
        int ind = 2*i;
        Real ds = wake.s[i+1] - wake.s[i];
        sigma_mWake[ind] = -1.0 / ds;
        sigma_mWake[ind+1] = 1.0 / ds;
    }
    
    // Ncoord value

    // Ncoords column
    for (int row=0; row<Ncoords+Nwake;++row){
        vsol.ue_m[colMajorIndex(row,Ncoords,Ncoords+Nwake)] = vsol.ue_sigma[colMajorIndex(row,Ncoords-1,Ncoords+Nwake)]*sigma_mWake[0];
    }


    // columns Ncoords+1 to end-1
    ue_mInt = (Ncoords+Nwake)*(Ncoords) + Ncoords+Nwake;
    sigma_mIndex = 1;
    for (int col=1; col<=Nwake-2;++col){

        Real sigma_mValue = sigma_mWake[sigma_mIndex];

        for (int row=0; row<Ncoords+Nwake; ++row){
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,Ncoords+col-2,Ncoords+Nwake)]*sigma_mValue;
        }

        sigma_mValue = sigma_mWake[sigma_mIndex+1];
        for (int row=0; row<Ncoords+Nwake; ++row){
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,Ncoords+col-1,Ncoords+Nwake)]*sigma_mValue;
        }

        sigma_mIndex += 2;
        ue_mInt += (Ncoords+Nwake);
    }

    // Last coluumn
    for (int row=0; row<Ncoords+Nwake;++row){
        vsol.ue_m[(Ncoords+Nwake)*(Ncoords+Nwake-1)+row] = vsol.ue_sigma[colMajorIndex(row,Ncoords+Nwake-3,Ncoords+Nwake)]*sigma_mWake[2*(Nwake-1)-1];
    }


    // sgnue switching
    for (int row = 0; row < Ncoords; ++row) {
        if (isol.edgeVelSign[row] == 1.0) {
            continue; // skip row — multiplying by 1 does nothing
        }
    
        for (int col = 0; col < Ncoords+Nwake; ++col) {
            vsol.ue_m[colMajorIndex(row,col,Ncoords+Nwake)] *= -1.0;
        }
    }

}