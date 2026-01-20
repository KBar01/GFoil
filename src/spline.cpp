#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.h"
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;
using namespace std::chrono;



struct CubicSpline1DOrig {
    int N = Nin;
    Real x[Nin] = {0};       // nodes (1D, N)
    Real c[4 * (Nin - 1)]={0};       // coefficients: (N-1) x 4 matrix stored column-major (4 rows)
};

struct CubicSpline1DFine {
    int N = Nfine;
    Real x[Nfine] = {0};       // nodes (1D, N)
    Real c[4 * (Nfine - 1)]={0};       // coefficients: (N-1) x 4 matrix stored column-major (4 rows)
};


struct Spline2DOrig {
    CubicSpline1DOrig Xspline;
    CubicSpline1DOrig Yspline;
};
struct Spline2DFine {
    CubicSpline1DFine Xspline;
    CubicSpline1DFine Yspline;
};


Real norm2(Real x, Real y) {
    return std::sqrt(x * x + y * y);
}

void fit_cubic_splineOrig(const Real* x, const Real* y, CubicSpline1DOrig& spline) {
    // Solve for natural cubic spline coefficients
    //spline.N = N;
    //spline.x = new Real[N];
    //spline.c = new Real[4 * (N - 1)];

    Real h[Nin - 1]={0};
    Real alpha[Nin - 1]={0};

    for (int i = 0; i < Nin; ++i){ spline.x[i] = x[i];}
    for (int i = 0; i < Nin - 1; ++i){
        h[i] = x[i+1] - x[i];
    }

    for (int i = 1; i < Nin - 1; ++i){
        alpha[i] = (3.0 / h[i]) * (y[i+1] - y[i]) - (3.0 / h[i-1]) * (y[i] - y[i-1]);
    }

    // Solve tridiagonal system for c coefficients
    Real l[Nin]={0}, mu[Nin]={0}, z[Nin]={0};
    l[0] = 1; mu[0] = 0; z[0] = 0;

    for (int i = 1; i < Nin - 1; ++i) {
        l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
    }

    l[Nin-1] = 1; z[Nin-1] = 0;

    Real c[Nin] = {0}; // temporary c coefficients
    Real b[Nin - 1] = {0};
    Real d[Nin - 1] = {0};
    

    for (int j = Nin - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j+1];
        b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j+1] - c[j]) / (3.0 * h[j]);

        // Save coefficients a=d[j], b=c[j], c=b[j], d=y[j]
        spline.c[IDX(0,j,4)] = d[j];
        spline.c[IDX(1,j,4)] = c[j];
        spline.c[IDX(2,j,4)] = b[j];
        spline.c[IDX(3,j,4)] = y[j];
    }
}


void fit_cubic_splineFine(const Real* x, const Real* y, CubicSpline1DFine& spline) {
    // Solve for natural cubic spline coefficients
    //spline.N = N;
    //spline.x = new Real[N];
    //spline.c = new Real[4 * (N - 1)];

    Real h[Nfine - 1]={0};
    Real alpha[Nfine - 1]={0};

    for (int i = 0; i < Nfine; ++i){ spline.x[i] = x[i];}
    for (int i = 0; i < Nfine - 1; ++i){
        h[i] = x[i+1] - x[i];
    }

    for (int i = 1; i < Nfine - 1; ++i){
        alpha[i] = (3.0 / h[i]) * (y[i+1] - y[i]) - (3.0 / h[i-1]) * (y[i] - y[i-1]);
    }

    // Solve tridiagonal system for c coefficients
    Real l[Nfine]={0}, mu[Nfine]={0}, z[Nfine]={0};
    l[0] = 1; mu[0] = 0; z[0] = 0;

    for (int i = 1; i < Nfine - 1; ++i) {
        l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
    }

    l[Nfine-1] = 1; z[Nfine-1] = 0;

    Real c[Nfine] = {0}; // temporary c coefficients
    Real b[Nfine - 1] = {0};
    Real d[Nfine - 1] = {0};

    for (int j = Nfine - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j+1];
        b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j+1] - c[j]) / (3.0 * h[j]);

        // Save coefficients a=d[j], b=c[j], c=b[j], d=y[j]
        spline.c[IDX(0,j,4)] = d[j];
        spline.c[IDX(1,j,4)] = c[j];
        spline.c[IDX(2,j,4)] = b[j];
        spline.c[IDX(3,j,4)] = y[j];
    }
}

void spline2dOrig(const Real* X, Spline2DOrig& spline) {
    
    Real xq[5] = {
        0.046910077030668, 0.230765344947158, 0.5,
        0.769234655052842, 0.953089922969332
    };
    Real wq[5] = {
        0.118463442528095, 0.239314335249683, 0.284444444444444,
        0.239314335249683, 0.118463442528095
    };

    // Allocate arc length arrays
    Real S[Nin]={0};
    Real Snew[Nin] = {0};

    // 1. Initial arc length estimate
    for (int i = 1; i < Nin; ++i) {
        Real dx = X[IDX(0,i,2)] - X[IDX(0,i-1,2)];
        Real dy = X[IDX(1,i,2)] - X[IDX(1,i-1,2)];
        S[i] = S[i-1] + std::sqrt(dx*dx + dy*dy);
    }

    // 2. Initial spline fit

    Real xPoints[Nin] = {0} ;
    Real yPoints[Nin] = {0};

    for (int i=0;i<Nin;++i){
        xPoints[i] = X[IDX(0,i,2)];
        yPoints[i] = X[IDX(1,i,2)];
    }

    fit_cubic_splineOrig(S, xPoints,spline.Xspline);
    fit_cubic_splineOrig(S, yPoints,spline.Yspline);

    // 3. Iterative correction
    constexpr int max_pass = 10;
    for (int pass = 0; pass < max_pass; ++pass) {
        Real serr = 0.0;
        Snew[0] = 0.0;

        for (int i = 0; i < Nin-1; ++i) {
            Real ds = S[i+1] - S[i];
            Real sint = 0.0;

            for (int q = 0; q < 5; ++q) {
                Real st = xq[q] * ds;
                Real xds = 3.0 * spline.Xspline.c[IDX(0,i,4)] * st * st +
                             2.0 * spline.Xspline.c[IDX(1,i,4)] * st +
                             spline.Xspline.c[IDX(2,i,4)];

                Real yds = 3.0 * spline.Yspline.c[IDX(0,i,4)] * st * st +
                             2.0 * spline.Yspline.c[IDX(1,i,4)] * st +
                             spline.Yspline.c[IDX(2,i,4)];

                sint += wq[q] * std::sqrt(xds*xds + yds*yds);
            }

            sint *= ds;
            serr = std::max(serr, std::abs(sint - ds));
            Snew[i+1] = Snew[i] + sint;
        }

        // Replace S with Snew
        for (int i = 0; i < Nin; ++i){
            S[i] = Snew[i];
        }

        // Refit splines
        fit_cubic_splineOrig(S, xPoints,spline.Xspline);
        fit_cubic_splineOrig(S, yPoints,spline.Yspline);
    }
}

void spline2dFine(const Real* X,Spline2DFine& spline) {
    
    const Real xq[5] = {
        0.046910077030668, 0.230765344947158, 0.5,
        0.769234655052842, 0.953089922969332
    };
    const Real wq[5] = {
        0.118463442528095, 0.239314335249683, 0.284444444444444,
        0.239314335249683, 0.118463442528095
    };

    // Allocate arc length arrays
    Real S[Nfine]={0};
    Real Snew[Nfine]={0};

    Real xPoints[Nfine]={0},yPoints[Nfine]={0};
    for (int i=0;i<Nfine;++i){
        xPoints[i] = X[IDX(0,i,2)];
        yPoints[i] = X[IDX(1,i,2)];
    }

    // 1. Initial arc length estimate
    for (int i = 1; i < Nfine; ++i) {
        Real dx = X[IDX(0,i,2)] - X[IDX(0,i-1,2)];
        Real dy = X[IDX(1,i,2)] - X[IDX(1,i-1,2)];
        S[i] = S[i-1] + std::sqrt(dx*dx + dy*dy);
    }

    // 2. Initial spline fit
    fit_cubic_splineFine(S, xPoints,spline.Xspline);
    fit_cubic_splineFine(S, yPoints,spline.Yspline);

    // 3. Iterative correction
    const int max_pass = 10;
    for (int pass = 0; pass < max_pass; ++pass) {
        Real serr = 0.0;
        Snew[0] = 0.0;

        for (int i = 0; i < Nfine-1; ++i) {
            Real ds = S[i+1] - S[i];
            Real sint = 0.0;

            for (int q = 0; q < 5; ++q) {
                Real st = xq[q] * ds;
                Real xds = 3.0 * spline.Xspline.c[IDX(0,i,4)] * st * st +
                             2.0 * spline.Xspline.c[IDX(1,i,4)] * st +
                             spline.Xspline.c[IDX(2,i,4)];

                Real yds = 3.0 * spline.Yspline.c[IDX(0,i,4)] * st * st +
                             2.0 * spline.Yspline.c[IDX(1,i,4)] * st +
                             spline.Yspline.c[IDX(2,i,4)];

                sint += wq[q] * std::sqrt(xds*xds + yds*yds);
            }

            sint *= ds;
            serr = std::max(serr, std::abs(sint - ds));
            Snew[i+1] = Snew[i] + sint;
        }

        // Replace S with Snew
        for (int i = 0; i < Nfine; ++i)
            S[i] = Snew[i];

        // Refit splines
        fit_cubic_splineFine(S, xPoints,spline.Xspline);
        fit_cubic_splineFine(S, yPoints,spline.Yspline);
    }
}

const Real xq[5] = {
    0.046910077030668, 0.230765344947158, 0.5,
    0.769234655052842, 0.953089922969332
};
const Real wq[5] = {
    0.118463442528095, 0.239314335249683, 0.284444444444444,
    0.239314335249683, 0.118463442528095
};

void evaluate_splineOrigtoFine(const CubicSpline1DOrig& spline, const Real* s, Real* yout) {
    
    for (int i = 0; i < Nfine; ++i) {
        Real val = s[i];
        int j = spline.N - 2;
        for (int k = 0; k < spline.N - 1; ++k)
            if (val < spline.x[k+1]) { j = k; break; }

        Real dx = val - spline.x[j];
        Real a = spline.c[IDX(0,j,4)];
        Real b = spline.c[IDX(1,j,4)];
        Real c = spline.c[IDX(2,j,4)];
        Real d = spline.c[IDX(3,j,4)];

        yout[i] = ((a * dx + b) * dx + c) * dx + d;
    }
}

void evaluate_splineFinetoOut(const CubicSpline1DFine& spline, const Real* s, Real* yout) {
    
    for (int i = 0; i < Ncoords; ++i) {
        Real val = s[i];
        int j = spline.N - 2;
        for (int k = 0; k < spline.N - 1; ++k)
            if (val < spline.x[k+1]) { j = k; break; }

        Real dx = val - spline.x[j];
        Real a = spline.c[IDX(0,j,4)];
        Real b = spline.c[IDX(1,j,4)];
        Real c = spline.c[IDX(2,j,4)];
        Real d = spline.c[IDX(3,j,4)];

        yout[i] = ((a * dx + b) * dx + c) * dx + d;
    }
}

Real cubic_interp1d(const Real* x, const Real* y, Real xq) {
    // Handle out-of-bounds queries
    if (xq <= x[0]) return y[0];
    if (xq >= x[Nfine - 1]) return y[Nfine - 1];

    // Find interval
    int i = 0;
    while (i < Nfine - 2 && xq > x[i + 1]) ++i;

    // Use 4 points: [i-1, i, i+1, i+2]
    int i0 = std::max(i - 1, 0);
    int i1 = i;
    int i2 = std::min(i + 1, Nfine - 1);
    int i3 = std::min(i + 2, Nfine - 1);

    // Normalize xq to [0, 1]
    Real t = (xq - x[i1]) / (x[i2] - x[i1]);

    // Cubic Hermite basis
    Real y0 = y[i0], y1 = y[i1], y2 = y[i2], y3 = y[i3];
    Real a0 = -0.5 * y0 + 1.5 * y1 - 1.5 * y2 + 0.5 * y3;
    Real a1 = y0 - 2.5 * y1 + 2.0 * y2 - 0.5 * y3;
    Real a2 = -0.5 * y0 + 0.5 * y2;
    Real a3 = y1;

    return ((a0 * t + a1) * t + a2) * t + a3;
}



void spline_curvature(
    const Real* Xin,    // Input points 2xNin
    Real Ufac, Real TEfac,     // Uniformity & TE resolution factors
    Real* Xout, Real (&Sout)[Ncoords]     // Output points (2xNcoords) and s values (1xNcoords)
) {
    // 1. Find X min/max
    Real xmin = Xin[0], xmax = Xin[0];
    for (int i = 1; i < Nin; ++i) {
        xmin = std::min(xmin, Xin[IDX(0,i,2)]);
        xmax = std::max(xmax, Xin[IDX(0,i,2)]);
    }

    // 2. Fit spline to original data
    Spline2DOrig PP;
    spline2dOrig(Xin, PP);

    // 3. Evaluate fine grid
    Real s[Nfine]={0};
    Real smax = PP.Xspline.x[PP.Xspline.N - 1];

    Real counter = 0 ;
    for (int i = 0; i < Nfine; ++i){
        s[i] = smax * counter / (Nfine - 1);
        counter += 1;
    }

    
    Real xfine[Nfine]={0};
    Real yfine[Nfine]={0};
    evaluate_splineOrigtoFine(PP.Xspline, s, xfine);
    evaluate_splineOrigtoFine(PP.Yspline, s, yfine);
    Real xyfine[2*Nfine]={0};
    for (int i=0;i<Nfine;++i){
        xyfine[IDX(0,i,2)] = xfine[i];
        xyfine[IDX(1,i,2)] = yfine[i];
    }

    // 4. Re-spline fine curve
    Spline2DFine PPfine;
    spline2dFine(xyfine, PPfine);

    for (int i=0;i<Nfine;++i){
        s[i] = PPfine.Xspline.x[i];
    }

    // 5. Compute curvature weights
    Real sk[Nfine]={0};

    for (int i = 0; i < Nfine-1; ++i) {
        Real ds = s[i+1] - s[i];
        Real sint = 0;

        for (int q = 0; q < 5; ++q) {
            Real st = s[i] + xq[q] * ds;
            int seg = PPfine.Xspline.N - 2;
            for (int j = 0; j < PPfine.Xspline.N - 1; ++j)
                if (st < PPfine.Xspline.x[j+1]) { seg = j; break; }

            Real dx2 = 6.0 * PPfine.Xspline.c[IDX(0,seg,4)] * (st - PPfine.Xspline.x[seg])
                       + 2.0 * PPfine.Xspline.c[IDX(1,seg,4)];
            Real dy2 = 6.0 * PPfine.Yspline.c[IDX(0,seg,4)] * (st - PPfine.Yspline.x[seg])
                       + 2.0 * PPfine.Yspline.c[IDX(1,seg,4)];

            sint += wq[q] * std::sqrt(dx2*dx2 + dy2*dy2);
        }

        sint = sint * ds * 0.5 + 0.01 * Ufac;

        Real xx = (0.5 * (xyfine[IDX(0,i,2)] + xyfine[IDX(0,i+1,2)]) - xmin) / (xmax - xmin);
        sint += TEfac * 0.5 * std::exp(-100 * (1.0 - xx));
        sk[i+1] = sk[i] + sint;
    }

    // 6. Avoid zero spacing
    Real sksum = 0 ;
    for (int i=0;i<Nfine;++i){sksum += sk[i];}
    
    for (int i = 0; i < Nfine; ++i){sk[i] += 2.0 * sksum / 501;}

    // 7. Map new s values
    Real skl[Ncoords]={0};
    for (int i = 0; i < Ncoords; ++i){
        skl[i] = sk[0] + i * (sk[Nfine - 1] - sk[0]) / (Ncoords - 1);
    }

    for (int i = 0; i < Ncoords; ++i) {
        Sout[i] = cubic_interp1d(sk, s, skl[i]);
    }

    // 8. Evaluate spline at new s values

    Real xXout[Ncoords]={0};
    Real yXout[Ncoords]={0};
    evaluate_splineFinetoOut(PPfine.Xspline, Sout, xXout);
    evaluate_splineFinetoOut(PPfine.Yspline, Sout, yXout);


    for (int i=0;i<Ncoords;++i){
        Xout[IDX(0,i,2)] = xXout[i];
        Xout[IDX(1,i,2)] = yXout[i];
    }
}

void make_panels(const Real (&inCoords)[2*Nin], Real (&outCoords)[2*Ncoords],const Real Ufac,const Real TEfac){


    Real outArcs[Ncoords];
    spline_curvature(inCoords,Ufac,TEfac,outCoords,outArcs);
};