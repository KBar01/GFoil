#pragma once

#include <codi.hpp>

inline int colMajorIndex(int row, int col, int num_rows) {
    return row + col*num_rows;
}

template<typename Real>
Real norm2(const Real* x) {
    return std::sqrt(x[0]*x[0] + x[1]*x[1]);
}
#define IDX(i,j,nrow) ((i)+(j)*(nrow)) // For col-major access
#define Nwake 30
#define Ncoords 200
constexpr int RVdimension = 4 * (Ncoords + Nwake);
#define Nfine 501
#define Nin 301
#define Nsound 250
#define NblPoints 250

