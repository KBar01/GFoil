#pragma once

// Shared template struct definitions used by both GFoil_fwd_codi and GFoil_AD.
// Nine structs are field-for-field identical between data_structs.h (fwd) and
// data_structs.hpp (AD); they live here as template<typename Real> with a _t
// suffix to avoid clashing with the using/typedef aliases in each build's own
// data_structs header.
//
// Caller must include real_type.h or real_type.hpp (for Ncoords, Nwake, norm2)
// and <vector> before including this header.

#include <cmath>
#include <vector>

template<typename Real>
struct Geom_t {
    Real chord = 1.0;
    Real wakelen = 1.0;
    int npoint = 1;
    Real xref[2] = {0.25, 0};
};

template<typename Real>
struct TE_t {
    Real t[2];
    Real hTE  = 0;
    Real dtdx = 0;
    Real tcp  = 0;
    Real tdp  = 0;

    TE_t(const Real* coords) { compute_TE_info(coords); }

    void compute_TE_info(const Real* coords) {
        Real lowerTang[2] = {coords[0] - coords[2], coords[1] - coords[3]};
        Real normLowerTang = norm2(lowerTang);
        lowerTang[0] /= normLowerTang;
        lowerTang[1] /= normLowerTang;

        Real upperTang[2] = {coords[2*(Ncoords-1)]   - coords[2*(Ncoords-2)],
                             coords[2*(Ncoords-1)+1]  - coords[2*(Ncoords-2)+1]};
        Real normUpperTang = norm2(upperTang);
        upperTang[0] /= normUpperTang;
        upperTang[1] /= normUpperTang;

        t[0] = 0.5*(lowerTang[0] + upperTang[0]);
        t[1] = 0.5*(lowerTang[1] + upperTang[1]);
        Real t_norm = norm2(t);
        t[0] /= t_norm;
        t[1] /= t_norm;

        Real TEVector[2] = {coords[2*(Ncoords-1)] - coords[0],
                            coords[2*(Ncoords-1)+1] - coords[1]};

        hTE  = -TEVector[0]*t[1] + TEVector[1]*t[0];
        dtdx = lowerTang[0]*upperTang[1] - upperTang[0]*lowerTang[1];

        Real norm_s = norm2(TEVector);
        Real p[2] = {TEVector[0]/norm_s, TEVector[1]/norm_s};

        tcp = std::abs(t[0]*p[1] - t[1]*p[0]);
        tdp = t[0]*p[0] + t[1]*p[1];
    }
};

template<typename Real>
struct Foil_t {
    Real x[2*Ncoords] = {0};
    Real s[Ncoords]   = {0};
    TE_t<Real> te;

    Foil_t(const Real* coordinates) : te(coordinates) {
        for (int i = 0; i < 2*Ncoords; i++)
            x[i] = coordinates[i];
        compute_arclengths();
    }

    void compute_arclengths() {
        s[0] = 0.0;
        for (int i = 1; i < Ncoords; ++i) {
            Real dx = x[2*i]   - x[2*(i-1)];
            Real dy = x[2*i+1] - x[2*(i-1)+1];
            s[i] = s[i-1] + std::sqrt(dx*dx + dy*dy);
        }
    }
};

template<typename Real>
struct Wake_t {
    Real x[Nwake*2] = {0};
    Real s[Nwake]   = {0};
    Real t[Nwake*2] = {0};
};

template<typename Real>
struct Oper_t {
    Real Vinf  = 1.0;
    Real alpha = 0.0;
    Real rho   = 1.225;
    Real Re    = 1e5;
    Real Ma    = 0.0;

    Oper_t(Real alpha_in, Real Re_in, Real Ma_in)
        : alpha(alpha_in), Re(Re_in), Ma(Ma_in) {}
};

template<typename Real>
struct Vsol_t {
    Real wgap[Nwake] = {0.0};
    Real ue_m[(Ncoords+Nwake)*(Ncoords+Nwake)] = {0.0};
    Real ue_sigma[(Ncoords+Nwake)*(Ncoords+Nwake-2)] = {0};
    bool turb[Ncoords+Nwake] = {false};
    std::vector<std::vector<int>> Is;
};

template<typename Real>
struct Post_t {
    Real cp[Ncoords+Nwake] = {0};
    Real cl = 0.0;
    Real cd = 0.0;
    Real cm = 0.0;
};

template<typename Real>
struct Param_t {
    Real rtol  = 1e-10;
    int  niglob = 50;

    Real ncrit = 9.0;
    Real Cuq   = 1.0;
    Real Dlr   = 0.9;
    Real SlagK = 5.6;

    Real CtauC = 1.8;
    Real CtauE = 3.3;

    Real GA = 6.7, GB = 0.75, GC = 18.0;

    Real Minf = 0.0, Vinf = 0.0, muinf = 0.0, mu0 = 0.0;
    Real rho0 = 1.0, H0   = 0.0, Tsrat = 0.35, gam = 1.4;
    Real KTb  = 1.0, KTl  = 0.0, cps   = 0.0;

    int breakLoop = 100;
};

template<typename Real>
struct Trans_t {
    int  isForced[2]  = {0, 0};
    int  transNode[2] = {1, 198};
    Real transPos[2]  = {0.5, 0.5};
};
