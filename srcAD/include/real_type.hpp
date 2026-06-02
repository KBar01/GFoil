#pragma once

#include <codi.hpp>

int colMajorIndex(int row, int col, int num_rows) {
    return row + col*num_rows;
}

template<typename Real>
Real norm2(const Real* x) {
    return std::sqrt(x[0]*x[0] + x[1]*x[1]);
}
#define IDX(i,j,nrow) ((i)+(j)*(nrow)) // For col-major access
#define Nwake 30
#define RVdimension 920
#define Ncoords 200
#define Nfine 501
#define Nin 301
#define Nsound 250
#define NblPoints 250

// Upper bound on Jacobian non-zeros: empirically ~14% of RVdimension² (≈130 entries
// per row of the 920×920 BL Jacobian). Observed peak NNZ is well below this limit.
static constexpr int RV_MAX_NNZ = 119700;

