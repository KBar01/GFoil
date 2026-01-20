
#include <codi.hpp>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "linsolveInv.hpp"


template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
 

int colMajorIndex(int row, int col, int num_rows) {
    return row + col*num_rows;
}

#ifndef USE_CODIPACK
    void solve_sys(Isol &isol, const Real* RHS) {
        int Nsize = Ncoords + 1;  // Matrix dimensions

        // Wrap the flattened column-major AIC vector into an Eigen matrix
        Eigen::Map<const Eigen::Matrix <Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> 
            A_eigen(isol.infMatrix, Nsize, Nsize);

        // Wrap the flattened column-major RHS vector into an Eigen matrix
        // Eigen-map is class that maps raw data (like c++ array) to eigen matrix/vector
        // WITHOUT copying data.
        // Inputs to this class are datatype (double in this case), dynamic means non-fixed at compile time
        // TODO: eigen sizes should be static.
        // ColMajor specifies how data is stored, matching arrays

        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> 
            rhs_eigen(RHS, Nsize, 2);
        
        // Solve Ax=b using Eigen's solver

        // Note: x here is a Eigen matrix 
        Eigen::MatrixXd x = A_eigen.colPivHouseholderQr().solve(rhs_eigen);

        // Ensure gamref is correctly mapped

        // This is taking the array of isol struct and temporarily mapping an
        // Eigen matrix to it, which allows easy update using solution as done in 
        // next step
        Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> 
            gamref_eigen(isol.gammasRef, Ncoords, 2);

        // Store solution in gamref
        gamref_eigen = x.block(0,0,Ncoords,2);

    }
#else
    void solve_sys(Isol &isol, const Real* RHS) {
        solve_sys_inv(isol,RHS);
    }
#endif


void build_gamma_codi(Isol &isol, const Foil& foil, const Oper& op) {
    // Build and solve the inviscid linear system for alpha=0,90,input
    // INPUT
    //   isol: Invisid solution structure
    //   foil: Aerofoil data  structure
    //   alpha: angle of attack (Radians)
    // OUTPUT
    //   isol.gamref: 0,90deg vorticity distributions at each node (Nx2)
    //   isol.gam: gamma for the particular input angle, alpha
    //   isol.infMatrix: aerodynamic influence coefficient matrix, filled in

    // Initialize matrices A and rhs
    Real rhs[2*(Ncoords + 1)] ; 

    // Initialise the panelinfo struct and influence coeffs used inside loop
    PanelInfo panelInfo;
    Real aij = 1.0;
    Real bij = 1.0;
    Real a   = 1.0;
    Real a_vortex = 1.0;
    Real b_vortex = 1.0;

    // Build influence matrix and rhs
    for (int i = 0; i < Ncoords; i++) {  // Loop over nodes

        for (int j = 0; j < Ncoords-1; j++) {  // Loop over panels
            // panel_linvortex_stream should be implemented to return aij and bij
            panel_linvortex_stream(foil.x[colMajorIndex(0,j,2)], foil.x[colMajorIndex(1,j,2)],
            foil.x[colMajorIndex(0,j+1,2)], foil.x[colMajorIndex(1,j+1,2)],
            foil.x[colMajorIndex(0,i,2)], foil.x[colMajorIndex(1,i,2)],
            panelInfo, aij, bij);
            
            isol.infMatrix[colMajorIndex(i,j,Ncoords+1)]   += aij;
            isol.infMatrix[colMajorIndex(i,j+1,Ncoords+1)] += bij;
        }

        isol.infMatrix[colMajorIndex(i,Ncoords,Ncoords+1)] = -1.0;  // Last unknown = streamfunction value on surface
        
        // Right-hand sides
        rhs[colMajorIndex(i,0,Ncoords+1)] = -foil.x[colMajorIndex(1,i,2)];
        rhs[colMajorIndex(i,1,Ncoords+1)] = foil.x[colMajorIndex(0,i,2)];

        // TE source influence
        panel_constsource_stream(foil.x[colMajorIndex(0,Ncoords-1,2)], foil.x[colMajorIndex(1,Ncoords-1,2)],
                              foil.x[colMajorIndex(0,0,2)], foil.x[colMajorIndex(1,0,2)],
                              foil.x[colMajorIndex(0,i,2)], foil.x[colMajorIndex(1,i,2)],
                              panelInfo, a);
        
        isol.infMatrix[colMajorIndex(i,0,Ncoords+1)] += -a*(0.5 * foil.te.tcp);
        isol.infMatrix[colMajorIndex(i,Ncoords-1,Ncoords+1)] += a*(0.5 * foil.te.tcp);

        // TE vortex panel 
        panel_linvortex_stream(foil.x[colMajorIndex(0,Ncoords-1,2)], foil.x[colMajorIndex(1,Ncoords-1,2)],
        foil.x[colMajorIndex(0,0,2)], foil.x[colMajorIndex(1,0,2)],
        foil.x[colMajorIndex(0,i,2)], foil.x[colMajorIndex(1,i,2)],
        panelInfo, a_vortex, b_vortex);
        
        isol.infMatrix[colMajorIndex(i,0,Ncoords+1)] += -(a_vortex + b_vortex)*(-0.5 * foil.te.tdp);
        isol.infMatrix[colMajorIndex(i,Ncoords-1,Ncoords+1)] += (a_vortex + b_vortex)*(-0.5 * foil.te.tdp);
    }

    // Kutta condition
    isol.infMatrix[colMajorIndex(Ncoords, 0, Ncoords+1)] = 1;
    isol.infMatrix[colMajorIndex(Ncoords, Ncoords-1, Ncoords+1)] = 1;
    //solve_sys(isol,rhs);   // Solves the linear sys using Eigen, and puts gammas in isol struct

    
    solve_sys(isol,rhs);


    for (int i = 0; i < Ncoords; ++i) {
        isol.gammas[i] = isol.gammasRef[colMajorIndex(i,0,Ncoords)]*std::cos(op.alpha) + isol.gammasRef[colMajorIndex(i,1,Ncoords)]*std::sin(op.alpha);
    }
};
