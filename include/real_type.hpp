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

// Default solver dimensions — override at runtime via struct constructors.
// Phase 2.3: Ncoords and Nwake are now default values; structs store their
// own copies so different calls can use different resolutions.
constexpr int Nwake   = 30;
constexpr int Ncoords = 200;
constexpr int RVdimension = 4 * (Ncoords + Nwake);  // legacy constant for AD fixed-size arrays
constexpr int Nfine   = 501;
constexpr int Nin     = 301;
constexpr int Nsound  = 250;
constexpr int NblPoints = 250;

