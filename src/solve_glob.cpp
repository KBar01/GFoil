#include <iostream>
#include <cmath>
#include <Eigen/Dense>


#include <Eigen/LU>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"


#include <vector>
#include <cmath>
#include <algorithm> // for std::max
#include <iostream>

#include <chrono>
#include <fstream>
#include "sparselinsolve.hpp"
using namespace std::chrono;




// === Drop-in sparse solver with full custom gradients for Codipack ===
// Place inside the #else branch where you currently have the dense QR solve.
// Requires Eigen (SparseLU) and Codipack headers already included.

#include <vector>
#include <iostream>
#include <cassert>

#include <vector>
#include <iostream>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>




void solve_sys(Glob& glob) {
    solve_sys_sparse(glob);
};


void solve_glob(const Foil&foil, const Isol&isol, Glob& glob, Vsol& vsol, const Oper& oper, const int doSolve) {
    
    
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
    get_ueinv(isol,ueinv);

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

    // all edge velocity indices
    int rowStart = 3*Nsys;
    for (int col=0;col<Nsys;++col){

        int colindex = 4*col + 3;
        for (int row = 0;row<Nsys;++row){
            
            //glob.R_V[colMajorIndex(rowStart+row,colindex,4*Nsys)] = (row == col ? 1.0 : 0.0) - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ds[col];
            Real zero = 0.0;
            glob.R_V_vals[glob.R_V_latest] = (row == col ? 1.0 : zero) - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ds[col];
            glob.R_V_rows[glob.R_V_latest] = rowStart+row;
            glob.R_V_cols[glob.R_V_latest] = colindex;
            glob.R_V_latest += 1 ;
        }
    }

    //all disp thickness indices
    for (int col=0;col<Nsys;++col){

        int colindex = 4*col + 1;
        for (int row = 0;row<Nsys;++row){
            //glob.R_V[colMajorIndex(rowStart+row,colindex,4*Nsys)] =  - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ue[col];
            glob.R_V_vals[glob.R_V_latest] = - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ue[col];
            glob.R_V_rows[glob.R_V_latest] = rowStart+row;
            glob.R_V_cols[glob.R_V_latest] = colindex;
            glob.R_V_latest += 1 ;
        }
    }

    if (doSolve) {
        solve_sys(glob);
    }
}


