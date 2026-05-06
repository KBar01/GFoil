#ifndef REAL_TYPE_H
#define REAL_TYPE_H



#include <codi.hpp>
using Real = codi::RealReverse;
using Tape = typename Real::Tape;

#define IDX(i,j,nrow) ((i)+(j)*(nrow)) // For col-major access
// Legacy header — forward solver uses unified include/real_type.hpp directly.
// These definitions kept for any remaining non-unified includes.
constexpr int Nwake   = 30;
constexpr int Ncoords = 200;
constexpr int RVdimension = 4 * (Ncoords + Nwake);
#define Nfine 501
#define Nin 301
#define Nsound 250
#define NblPoints 250

#endif