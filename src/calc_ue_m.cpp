#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <codi.hpp>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "vector_ops.hpp"
#include "linSolveMatrix.hpp"


void solve_sys(Isol&isol, const Real* RHS, Real*Bp){
    solve_sys_ue(isol,RHS,Bp);
}

void compute_Dw(const Real* Cgam, const Real* Bp, const Real* Csig, Real* Dw){

    constexpr int nPanels = Nwake+Ncoords-2 ;
    cnp::matmat_mul<Nwake,Ncoords,nPanels>(Cgam,Bp,Dw);
    cnp::add_inplace<(Nwake*nPanels)>(Dw,Csig);
}

void calc_ue_m(const Foil&foil,const Wake&wake,Isol&isol,Vsol &vsol) {


    // Allocate Cgam: Nwake x Ncoords
    Real Cgam[Nwake * Ncoords] = {0.0}; 
    for (int i = 0; i < Nwake; ++i) {
        
        Real V_G[2*Ncoords] = {0.0};
        dvelocity_dgamma(foil,wake.x[colMajorIndex(0,i,2)],wake.x[colMajorIndex(1,i,2)],V_G);

        for (int j = 0; j < Ncoords; ++j) {
            Real V_Gx = V_G[colMajorIndex(0, j, 2)];
            Real V_Gy = V_G[colMajorIndex(1, j, 2)];
            Real Cgamx = V_Gx*wake.t[colMajorIndex(0, i, 2)];
            Real Cgamy = V_Gy*wake.t[colMajorIndex(1, i, 2)];
            Cgam[colMajorIndex(i, j, Nwake)] = Cgamx+Cgamy;
        }
    }

    constexpr int nPanels = Ncoords + Nwake - 2;

    // Allocate B: (Ncoords+1) x npan
    Real B[(Ncoords+1) * nPanels] = {0.0};

    PanelInfo info;
    for (int i = 0; i < Ncoords; ++i) {
        
        Real foilPoint[2] = {foil.x[colMajorIndex(0, i, 2)], foil.x[colMajorIndex(1, i, 2)]};

        // Airfoil panels
        for (int j = 0; j < Ncoords-1; ++j) {
            
            Real foilPanelx1 = foil.x[colMajorIndex(0, j, 2)];
            Real foilPanely1 = foil.x[colMajorIndex(1, j, 2)];
            Real foilPanelx2 = foil.x[colMajorIndex(0, j+1, 2)];
            Real foilPanely2 = foil.x[colMajorIndex(1, j+1, 2)];
            Real a;
            panel_constsource_stream(foilPanelx1,foilPanely1,foilPanelx2,foilPanely2,foilPoint[0],foilPoint[1],info,a);
            B[colMajorIndex(i, j, Ncoords+1)] = a;
        }

        // Wake panels
        for (int j = 0; j < Nwake-1; ++j) {
            
            Real wakePanelLeft[2] = {wake.x[colMajorIndex(0,j,2)],wake.x[colMajorIndex(1,j,2)]};  // 3 pts: left, mid, right
            Real wakePanelRight[2] = {wake.x[colMajorIndex(0,j+1,2)],wake.x[colMajorIndex(1,j+1,2)]};  // 3 pts: left, mid, right
            Real wakePanelMid[2] = {0.5*(wakePanelLeft[0]+wakePanelRight[0]),0.5*(wakePanelLeft[1]+wakePanelRight[1])};
            
            if (j == Nwake-2) {
    
                wakePanelRight[0] = 2*wakePanelRight[0]-wakePanelMid[0];  // ghost extension
                wakePanelRight[1] = 2*wakePanelRight[1]-wakePanelMid[1];
            }
            
            Real a, b;

            // Left half
            panel_linsource_stream(wakePanelLeft[0],wakePanelLeft[1],wakePanelMid[0],wakePanelMid[1],foilPoint[0],foilPoint[1],info,a,b);
            if (j > 0) {
                B[colMajorIndex(i, Ncoords-1+j, Ncoords+1)] += 0.5*a + b;
                B[colMajorIndex(i, Ncoords-1+j-1, Ncoords+1)] += 0.5*a;
            } else {
                B[colMajorIndex(i, Ncoords-1+j, Ncoords+1)] += b;
            }

            // Right half
            panel_linsource_stream(wakePanelMid[0],wakePanelMid[1],wakePanelRight[0],wakePanelRight[1],foilPoint[0],foilPoint[1],info,a,b);
            B[colMajorIndex(i, Ncoords-1+j, Ncoords+1)] += a + 0.5*b;
            if (j < Nwake-2) {
                B[colMajorIndex(i, Ncoords+j, Ncoords+1)] += 0.5*b;
            } else {
                B[colMajorIndex(i, Ncoords-1+j, Ncoords+1)] += 0.5*b;
            }
        }
    }

    // Solve Bp = -AIC^{-1} * B (note B is (Ncoords+1)xnpan)
    Real Bp[Ncoords * nPanels];
    solve_sys(isol,B,Bp);

    // Csig: Nwake x npan
    Real Csig[Nwake * nPanels]={0.0};

    Real a1,a2,b1,b2;
    for (int i = 0; i < Nwake; ++i) {
        
        Real xi[2] = {wake.x[colMajorIndex(0, i, 2)],wake.x[colMajorIndex(1, i, 2)]};
        Real ti[2] = {wake.t[colMajorIndex(0, i, 2)],wake.t[colMajorIndex(1, i, 2)]};

        // Constant sources on foil
        int jstart = (i == 0) ? 1 : 0;
        int jend   = (i == 0) ? Ncoords-2 : Ncoords-1;

        for (int j = jstart; j < jend; ++j) {
            
            panel_constsource_velocity(
                foil.x[colMajorIndex(0,j,2)],
                foil.x[colMajorIndex(1,j,2)],
                foil.x[colMajorIndex(0,j+1,2)],
                foil.x[colMajorIndex(1,j+1,2)],
                xi[0],xi[1],info,a1,a2);
            
            Csig[colMajorIndex(i, j, Nwake)] = a1*ti[0] + a2*ti[1];

        }
        // piecewise linear sources across wake panel halves (else singular)
        for (int j = 0; j < Nwake; ++j) {
            
            
            int I0 = std::max(j-1, 0);
            int I1 = j;
            int I2 = std::min(j+1, Nwake-1);
        
            Real leftX = 0.5*(wake.x[colMajorIndex(0,I0,2)] + wake.x[colMajorIndex(0,I1,2)]);
            Real leftY = 0.5*(wake.x[colMajorIndex(1,I0,2)] + wake.x[colMajorIndex(1,I1,2)]);
            Real midX  =      wake.x[colMajorIndex(0,I1,2)];
            Real midY  =      wake.x[colMajorIndex(1,I1,2)];
            Real rightX = 0.5*(midX + wake.x[colMajorIndex(0,I2,2)]);
            Real rightY = 0.5*(midY + wake.x[colMajorIndex(1,I2,2)]);

        
            if (j == Nwake - 1) {
                rightX = 2*midX - leftX;
                rightY = 2*midY - leftY;
            }
            
            Real leftLength[2] = {
                midX - leftX,
                midY - leftY
            };

            Real rightLength[2] = {
                rightX - midX,
                rightY - midY
            };

            Real d1 = norm2(leftLength); // norm between center and mid-left
            Real d2 = norm2(rightLength); // norm between mid-right and center
        
            if (i == j) {
                // first point: special TE system (three panels meet)
                if (j == 0) {
                    
                    Real lowPanelLength[2] = {
                        foil.x[colMajorIndex(0,1,2)]-foil.x[colMajorIndex(0,0,2)],
                        foil.x[colMajorIndex(1,1,2)]-foil.x[colMajorIndex(1,0,2)]
                    };

                    Real upPanelLength[2] = {
                        foil.x[colMajorIndex(0,Ncoords-1,2)]-foil.x[colMajorIndex(0,Ncoords-2,2)],
                        foil.x[colMajorIndex(1,Ncoords-1,2)]-foil.x[colMajorIndex(1,Ncoords-2,2)]
                    };

                    Real dl = norm2(lowPanelLength);
                    Real du = norm2(upPanelLength);
        
                    Csig[colMajorIndex(i, 0, Nwake)]         += (0.5/M_PI) * (std::log(dl/d2) + 1.0);
                    Csig[colMajorIndex(i, Ncoords-2, Nwake)] += (0.5/M_PI) * (std::log(du/d2) + 1.0);
                    Csig[colMajorIndex(i, Ncoords-1, Nwake)] += -0.5/M_PI;
                
                } else if (j == Nwake - 1) {
                    // self effect = 0
                } else {
                    Real aa = (0.25/M_PI) * std::log(d1/d2);
                    Csig[colMajorIndex(i, Ncoords-2+j, Nwake)] += aa + 0.5/M_PI;
                    Csig[colMajorIndex(i, Ncoords-1+j, Nwake)] += aa - 0.5/M_PI;
                }
            } else {
                if (j == 0) {
                    
                    
                    panel_linsource_velocity(midX,midY,rightX,rightY,xi[0],xi[1],info,a1,b1,a2,b2);

                    Real a = a1*ti[0] + a2*ti[1];
                    Real b = b1*ti[0] + b2*ti[1];

                    Csig[colMajorIndex(i, Ncoords-1, Nwake)]   += b;
                    Csig[colMajorIndex(i, 0, Nwake)]           += a;
                    Csig[colMajorIndex(i, Ncoords-2, Nwake)]   += a;

                } else if (j == Nwake-1) {
                    
                    panel_constsource_velocity(leftX,leftY,rightX,rightY,xi[0],xi[1],info,a1,a2); // Xj[:, [0,2]]
                    Csig[colMajorIndex(i, Ncoords+Nwake-3, Nwake)] += (a1*ti[0] + a2*ti[1]);

                } else {

                    panel_linsource_velocity(leftX,leftY,midX,midY,xi[0],xi[1],info,a1,b1,a2,b2);
                    Csig[colMajorIndex(i, Ncoords-2+j, Nwake)]     += (a1*ti[0] + a2*ti[1]) + 0.5 * (b1*ti[0] + b2*ti[1]);
                    
                    panel_linsource_velocity(midX,midY,rightX,rightY,xi[0],xi[1],info,a1,b1,a2,b2);   // Xj[:, [1,2]]
                    Csig[colMajorIndex(i, Ncoords - 1 + j, Nwake)] += 0.5 * (a1*ti[0] + a2*ti[1]) + (b1*ti[0] + b2*ti[1]);
                }
            }
        }
    }

    // Combine Dw = Cgam * Bp + Csig
    Real Dw[Nwake * nPanels];
    //compute_Dw_eigen(Cgam,Bp,Csig,Dw);  // Cgam [Nwake x Ncoords], Bp [Ncoords x npan], Dw [Nwake x npan]

    compute_Dw(Cgam,Bp,Csig,Dw);
    for (int j = 0; j < nPanels; ++j) {
        Dw[colMajorIndex(0, j, Nwake)] = Bp[colMajorIndex(Ncoords-1, j, Ncoords)];
    }


    for (int i=0; i<Ncoords;++i){
        for (int j=0; j<nPanels; ++j){
            vsol.ue_sigma[colMajorIndex(i,j,Ncoords+Nwake)] = Bp[colMajorIndex(i,j,Ncoords)];
        }
    }

    for (int i=0; i<Nwake;++i){
        for (int j=0; j<nPanels; ++j){
            vsol.ue_sigma[colMajorIndex(Ncoords+i,j,Ncoords+Nwake)] = Dw[colMajorIndex(i,j,Nwake)];
        }
    }

}

