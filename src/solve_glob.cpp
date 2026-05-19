#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "real_type.h"
#include "panel_funcs.hpp"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"
#include "solver_funcs.hpp"

#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include "sparselinsolve.hpp"
using namespace std::chrono;

#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>


void solve_sys(Glob& glob) {
    solve_sys_sparse(glob);
};


void solve_glob(const Foil&foil, const Isol&isol, Glob& glob, Vsol& vsol, const Oper& oper, const int doSolve) {

    constexpr int Nsys = Ncoords+Nwake;

    // Steps 1-3: shared kernel fills glob.R[3*Nsys:]
    // Pass isol twice: it plays both the "isolc" (gammas/uewi) and
    // "isolv" (edgeVelSign) roles since fwd Isol is a unified struct.
    ue_residual_kernel<Real>(isol, isol, vsol, glob);

    // Fwd-only: sparse Jacobian fill for ue rows.
    // Recompute ue[] and ds[] from glob.U (same values the kernel used).
    Real uemax = 0.0;
    for (int i = 0; i < Nsys; ++i)
        uemax = std::max(uemax, std::abs(glob.U[colMajorIndex(3,i,4)]));

    Real ue[Nsys];
    Real ds[Nsys];
    for (int i = 0; i < Nsys; ++i) {
        ue[i] = std::max(glob.U[colMajorIndex(3,i,4)], Real(1e-10)*uemax);
        ds[i] = glob.U[colMajorIndex(1,i,4)];
    }

    int rowStart = 3*Nsys;

    // d(R_ue)/d(ue) columns
    for (int col=0;col<Nsys;++col){
        int colindex = 4*col + 3;
        for (int row = 0;row<Nsys;++row){
            Real zero = 0.0;
            glob.R_V_vals[glob.R_V_latest] = (row == col ? 1.0 : zero) - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ds[col];
            glob.R_V_rows[glob.R_V_latest] = rowStart+row;
            glob.R_V_cols[glob.R_V_latest] = colindex;
            glob.R_V_latest += 1;
        }
    }

    // d(R_ue)/d(ds) columns
    for (int col=0;col<Nsys;++col){
        int colindex = 4*col + 1;
        for (int row = 0;row<Nsys;++row){
            glob.R_V_vals[glob.R_V_latest] = -vsol.ue_m[colMajorIndex(row,col,Nsys)]*ue[col];
            glob.R_V_rows[glob.R_V_latest] = rowStart+row;
            glob.R_V_cols[glob.R_V_latest] = colindex;
            glob.R_V_latest += 1;
        }
    }

    if (doSolve) {
        solve_sys(glob);
    }
}
