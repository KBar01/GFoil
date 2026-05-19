#ifndef REAL_TYPE_H
#define REAL_TYPE_H



#include <codi.hpp>
using Real = codi::RealReverse;
using Tape = typename Real::Tape;

#include <cmath>
template<typename T> inline T norm2(const T* x)    { return std::sqrt(x[0]*x[0]+x[1]*x[1]); }
template<typename T> inline T norm2_3D(const T* x) { return std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }
inline int colMajorIndex(int row, int col, int num_rows) {
    return row + col * num_rows;
}

#define IDX(i,j,nrow) ((i)+(j)*(nrow)) // For col-major access
#define Nwake 30
#define RVdimension 920
#define Ncoords 200
#define Nfine 501
#define Nin 301
#define Nsound 250
#define NblPoints 250

#endif