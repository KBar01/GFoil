#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"

int march_amplification(Glob &glob, Vsol &vsol, Isol &isol, int si, const Param&param, Trans&tdata, const bool force) {
    
    
    const std::vector<int> &Is = vsol.Is[si];
    int N = Is.size();

    glob.U[colMajorIndex(2,Is[0],4)] = 0.0; // initial amplification

    int i = 1;
    while (i < N) {
        
        int i1 = Is[i-1];
        int i2 = Is[i];

        Real U1[4], U2[4];      // these are copies to avoid changing
        for (int j = 0; j < 4; ++j) {
            U1[j] = glob.U[colMajorIndex(j,i1,4)];
            U2[j] = glob.U[colMajorIndex(j,i2,4)];
        }
        if (vsol.turb[i2]) {U2[2] = U1[2]*1.01;}
        
        Real dx = isol.distFromStag[i2] - isol.distFromStag[i1];

        constexpr int nNewton = 20;

        const Real one = 1, zero=0 ;

        for (int iNewton=0; iNewton < nNewton; ++iNewton) {
             
            Real damp1, damp2, damp1_U[4], damp2_U[4];
            damp1 = get_damp(U1[0],U1[1],U1[2],U1[3],param,damp1_U);
            damp2 = get_damp(U2[0],U2[1],U2[2],U2[3],param,damp2_U);

            Real damp, damp_U[8];
            damp = upwind_half(damp1, damp1_U, damp2, damp2_U, damp_U);

            Real Ramp = U2[2] - U1[2] - damp*dx;

            if (std::abs(Ramp) < 1e-12) break;

            Real Ramp_U[8] = {0.0};
            Ramp_U[2] = -1.0;
            Ramp_U[6] = 1.0;
            for (int j = 0; j < 8; ++j){ Ramp_U[j] -= damp_U[j]*dx;}

            Real dU = -Ramp / Ramp_U[6];
            Real dmax = 0.5 * (1.01 - static_cast<Real>(iNewton)/nNewton);
            Real omega = (std::abs(dU) > dmax) ? dmax/std::abs(dU) : one;
            U2[2] += omega*dU;
        }
        
        // OR if node is a forced transition node break also

        if ((U2[2] > param.ncrit)) {
            tdata.isForced[si] = 0 ;
            break;
        }
        else if (force && (i1==tdata.transNode[si])){
            
            tdata.isForced[si] = 1 ;
            break;
        }
        else {
            glob.U[colMajorIndex(2,i2,4)] = U2[2];
        }
        ++i;
    }

    return i - 1;
}


void update_transition(Glob &glob, Vsol &vsol, Isol &isol, Param&param, Trans&tdata,const bool force) {
    
    
    for (int si = 0; si < 2; ++si) {
        

        const std::vector<int> &Is = vsol.Is[si];
        int nSurfPoints = Is.size();
        
        
        // find current last laminar station
        int ilam0 = nSurfPoints - 1;
        for (int i = 0; i < nSurfPoints; ++i) {
            if (vsol.turb[Is[i]]) {
                ilam0 = i - 1;
                break;
            }
        }

        // copy current amp/shear
        Real sa[Ncoords];
        for (int state=0;state<Ncoords;++state){
            sa[state] = glob.U[colMajorIndex(2,state,4)];
        }
    
        // get new transition location
        int ilam = march_amplification(glob, vsol, isol, si, param,tdata,force);

        if (ilam == ilam0) {
            
            for (int state=0;state<Ncoords;++state){
                glob.U[colMajorIndex(2,state,4)] = sa[state];
            }
            continue; // no change to states in this case
        }

        if (ilam < ilam0) {
            bool turb = true;
            Real sa0, cttr_U[4];
            sa0 = get_cttr(glob.U[colMajorIndex(0,Is[ilam+1],4)],
                glob.U[colMajorIndex(1,Is[ilam+1],4)],
                glob.U[colMajorIndex(2,Is[ilam+1],4)],
                glob.U[colMajorIndex(3,Is[ilam+1],4)],
                turb,param,cttr_U);

            Real sa1 = (ilam0 < nSurfPoints-1) ? glob.U[colMajorIndex(2,Is[ilam0+1],4)] : sa0;

            const Real zero=0,one=1;
        
            Real xi_start = isol.distFromStag[Is[ilam+1]];
            Real xi_end = isol.distFromStag[Is[std::min(ilam0+1, nSurfPoints-1)]];
            Real dx = xi_end - xi_start;

            for (int i = ilam+1; i <= ilam0; ++i) {
                
                Real f = (dx == 0 || i == ilam+1) ? zero : (isol.distFromStag[Is[i]]-xi_start)/dx;
                
                if ((ilam+1) == ilam0) f = one;

                Real sa_interp = sa0 + f*(sa1 - sa0);
    
                glob.U[colMajorIndex(2,Is[i],4)] = sa_interp;
                vsol.turb[Is[i]] = true;
            }
        }
        else if (ilam > ilam0){

            for (int i = ilam0; i <= ilam; ++i)
                vsol.turb[Is[i]] = false;
        }
 
    }
}

