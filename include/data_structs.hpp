#pragma once

#include <vector>
#include <string>
#include <cmath>
#include "real_type.hpp"

//-------------------------------------------------------------------------------
// Geometry Struct
template<typename Real>
struct Geom {
    Real chord   = 1.0;
    Real wakelen = 1.0;
    int  npoint  = 1;
    Real xref[2] = {0.25, 0};
};

//-------------------------------------------------------------------------------
// TE Struct
template<typename Real>
struct TE {
    Real t[2];
    Real hTE = 0, dtdx = 0, tcp = 0, tdp = 0;

    TE(const Real* coords, int nc) { compute_TE_info(coords, nc); }

    void compute_TE_info(const Real* coords, int nc) {
        Real lowerTang[2] = {coords[0] - coords[2], coords[1] - coords[3]};
        Real n1 = norm2(lowerTang); lowerTang[0]/=n1; lowerTang[1]/=n1;

        Real upperTang[2] = {coords[2*(nc-1)] - coords[2*(nc-2)],
                             coords[2*(nc-1)+1] - coords[2*(nc-2)+1]};
        Real n2 = norm2(upperTang); upperTang[0]/=n2; upperTang[1]/=n2;

        t[0] = 0.5*(lowerTang[0]+upperTang[0]);
        t[1] = 0.5*(lowerTang[1]+upperTang[1]);
        Real tn = norm2(t); t[0]/=tn; t[1]/=tn;

        Real TEV[2] = {coords[2*(nc-1)]-coords[0], coords[2*(nc-1)+1]-coords[1]};
        hTE  = -TEV[0]*t[1] + TEV[1]*t[0];
        dtdx = lowerTang[0]*upperTang[1] - upperTang[0]*lowerTang[1];
        Real ns = norm2(TEV);
        Real p[2] = {TEV[0]/ns, TEV[1]/ns};
        tcp = std::abs(t[0]*p[1] - t[1]*p[0]);
        tdp = t[0]*p[0] + t[1]*p[1];
    }
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Foil {
    int nc;                   // number of surface nodes
    std::vector<Real> x;      // flattened 2D coords (2*nc)
    std::vector<Real> s;      // arc-length at each node (nc)
    TE<Real> te;

    Foil(const Real* coordinates, int nc_ = Ncoords)
        : nc(nc_), x(2*nc_), s(nc_), te(coordinates, nc_)
    {
        for (int i=0; i<2*nc_; ++i) x[i] = coordinates[i];
        compute_arclengths();
    }

    void compute_arclengths() {
        s[0] = 0;
        for (int i=1; i<nc; ++i) {
            Real dx = x[2*i]   - x[2*(i-1)];
            Real dy = x[2*i+1] - x[2*(i-1)+1];
            s[i] = s[i-1] + std::sqrt(dx*dx + dy*dy);
        }
    }
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Wake {
    int nw;
    std::vector<Real> x;  // 2*nw
    std::vector<Real> s;  // nw
    std::vector<Real> t;  // 2*nw

    Wake(int nw_ = Nwake) : nw(nw_), x(2*nw_,0), s(nw_,0), t(2*nw_,0) {}
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Oper {
    Real Vinf = 1.0, alpha = 0.0, rho = 1.225, Re = 1e5, Ma = 0.0;
    Oper(Real a, Real Re_, Real Ma_) : alpha(a), Re(Re_), Ma(Ma_) {}
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Isolc {
    int nc, nw;
    std::vector<Real> infMatrix;  // (nc+1)^2
    std::vector<Real> gammasRef;  // 2*nc
    std::vector<Real> gammas;     // nc
    std::vector<Real> uewi;       // nw
    std::vector<Real> uewiref;    // 2*nw

    Isolc(int nc_ = Ncoords, int nw_ = Nwake)
        : nc(nc_), nw(nw_),
          infMatrix((nc_+1)*(nc_+1), 0),
          gammasRef(2*nc_, 0), gammas(nc_, 0),
          uewi(nw_, 0), uewiref(2*nw_, 0) {}
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Isolv {
    int nc, nw;
    Real stagArcLocation = 0;
    Real stagXLocation[2] = {0, 0};
    Real sstag_g[2]  = {0, 0};
    Real sstag_ue[2] = {0, 0};
    int  stagIndex[2] = {0, 0};
    std::vector<int>  edgeVelSign;  // nc
    std::vector<Real> distFromStag; // nc+nw

    Isolv(int nc_ = Ncoords, int nw_ = Nwake)
        : nc(nc_), nw(nw_),
          edgeVelSign(nc_, -1),
          distFromStag(nc_+nw_, 0) {}
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Vsol {
    int nc, nw;
    std::vector<Real> wgap;      // nw
    std::vector<Real> ue_sigma;  // (nc+nw)*(nc+nw-2)
    std::vector<Real> ue_m;      // (nc+nw)^2
    std::vector<int>  turb;      // nc+nw  (0=laminar, 1=turbulent; int avoids std::vector<bool>)
    std::vector<std::vector<int>> Is;

    Vsol(int nc_ = Ncoords, int nw_ = Nwake)
        : nc(nc_), nw(nw_),
          wgap(nw_, 0),
          ue_sigma((nc_+nw_)*(nc_+nw_-2), 0),
          ue_m((nc_+nw_)*(nc_+nw_), 0),
          turb(nc_+nw_, 0) {}
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Glob {
    int nc, nw;
    std::vector<Real> U;   // 4*(nc+nw)
    std::vector<Real> dU;  // 4*(nc+nw)
    std::vector<Real> R;   // 4*(nc+nw)

    std::vector<Real> R_V_vals;
    std::vector<int>  R_V_rows;
    std::vector<int>  R_V_cols;

    int  convergenceIteration = -1;
    bool doADrestartExtract   = false;
    Real ADrestartRnorm       = 1.0e-10;

    Glob(int nc_ = Ncoords, int nw_ = Nwake)
        : nc(nc_), nw(nw_),
          U(4*(nc_+nw_), 0),
          dU(4*(nc_+nw_), 0),
          R(4*(nc_+nw_), 0) {}
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Post {
    int nc, nw;
    std::vector<Real> cp;  // nc+nw
    Real cl = 0, cd = 0, cm = 0;

    Post(int nc_ = Ncoords, int nw_ = Nwake)
        : nc(nc_), nw(nw_), cp(nc_+nw_, 0) {}
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Param {
    Real rtol  = 1e-10;
    int  niglob = 50;
    Real ncrit  = 9.0, Cuq = 1.0, Dlr = 0.9, SlagK = 5.6;
    Real CtauC  = 1.8, CtauE = 3.3;
    Real GA = 6.7, GB = 0.75, GC = 18.0;
    Real Minf = 0, Vinf = 0, muinf = 0, mu0 = 0;
    Real rho0 = 1, H0 = 0, Tsrat = 0.35, gam = 1.4;
    Real KTb = 1, KTl = 0, cps = 0;
    int breakLoop = 100;
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Trans {
    int  isForced[2]  = {0, 0};
    int  transNode[2] = {1, 198};
    Real transPos[2]  = {0.5, 0.5};
};
