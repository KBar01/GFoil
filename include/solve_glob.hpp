#pragma once
#include "real_type.hpp"
#include "data_structs.hpp"
#include "get_funcs.hpp"
#ifdef USE_CODIPACK
#include "sparselinsolve.hpp"
#endif

template<typename Real>
void solve_glob(const Foil<Real>& foil, const Isolc<Real>& isolc, const Isolv<Real>& isol,
    Glob<Real>& glob, Vsol<Real>& vsol, const Oper<Real>& oper, const int doSolve)
{
    const int Nsys = glob.nc + glob.nw;

    // Step 1: Clamp edge velocity to avoid near-zero entries in the ds Jacobian block.
    // Near-zero ue values at the stagnation point make ue_m*ue -> 0, causing a
    // numerically singular matrix. The clamp matches the original src/solve_glob.cpp.
    std::vector<Real> ue(Nsys, 0), ds(Nsys, 0);
    Real uemax = 0;
    for (int i=0; i<Nsys; ++i)
        uemax = std::max(uemax, std::abs(glob.U[colMajorIndex(3,i,4)]));
    for (int i=0; i<Nsys; ++i) {
        ue[i] = std::max(glob.U[colMajorIndex(3,i,4)], Real(1e-10)*uemax);
        ds[i] = glob.U[colMajorIndex(1,i,4)];
    }

    // Step 2: Compute inviscid edge velocity
    std::vector<Real> ueinv(Nsys, 0);
    get_ueinv<Real>(isolc, isol, ueinv.data());

    // Step 3: Set edge-velocity residuals glob.R[3*Nsys..4*Nsys-1]
    //   R_ue[i] = ue[i] - ueinv[i] - (ue_m * (ds * ue))[i]
    std::vector<Real> tempRHS(Nsys, 0);
    for (int i=0; i<Nsys; ++i) tempRHS[i] = ds[i]*ue[i];

    Real* Rpointer = &glob.R[3*Nsys];
    // Rpointer[i] = sum_j ue_m[i,j] * tempRHS[j]
    for (int i=0; i<Nsys; ++i) {
        Real s = 0;
        for (int j=0; j<Nsys; ++j)
            s += vsol.ue_m[colMajorIndex(i,j,Nsys)] * tempRHS[j];
        Rpointer[i] = s;
    }
    for (int i=0; i<Nsys; ++i)
        Rpointer[i] = ue[i] - (ueinv[i] + Rpointer[i]);

    // Step 4: Assemble ue-equation Jacobian rows (3*Nsys .. 4*Nsys-1)
    const int rowStart = 3*Nsys;
    for (int col=0; col<Nsys; ++col) {
        const int colindex_ue = 4*col+3;
        const int colindex_ds = 4*col+1;
        for (int row=0; row<Nsys; ++row) {
            const Real ue_m_rc = vsol.ue_m[colMajorIndex(row,col,Nsys)];
            // d(R_ue_row)/d(ue_col)
            glob.R_V_vals.push_back((row==col ? Real(1.0) : Real(0.0)) - ue_m_rc*ds[col]);
            glob.R_V_rows.push_back(rowStart+row);
            glob.R_V_cols.push_back(colindex_ue);
            // d(R_ue_row)/d(ds_col)
            glob.R_V_vals.push_back(-ue_m_rc * ue[col]);
            glob.R_V_rows.push_back(rowStart+row);
            glob.R_V_cols.push_back(colindex_ds);
        }
    }

    if (doSolve) {
        // Diagnostic: check ue_m diagonal and ue values
        double ue_m_diag0 = vsol.ue_m[0].getValue();
        double ue0 = ue[0].getValue();
        double ds0 = ds[0].getValue();
        double R_ue0 = glob.R[3*Nsys].getValue();
        fprintf(stderr, "DBG solve_glob: ue_m[0,0]=%.4g ue[0]=%.4g ds[0]=%.4g R_ue[0]=%.4g nnz_before=%d\n",
            ue_m_diag0, ue0, ds0, R_ue0, (int)glob.R_V_vals.size());
        solve_sys_sparse<Real>(glob);
    }
}
