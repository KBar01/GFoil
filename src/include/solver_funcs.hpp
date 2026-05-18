#pragma once

// Shared template implementations used by both GFoil_fwd_codi and
// GFoil_AD. Replaces duplicated non-template functions that previously
// existed in src/*.cpp and srcAD/include/main_func.hpp.
// Duck-typed on struct params (FoilT, IsolT, etc.) so it works with
// both the non-template fwd structs and the template AD structs.

#include <cmath>
#include <vector>

// ── range (non-template helper) ──────────────────────────────────────────────
inline std::vector<int> range(int start, int end, int step = 1) {
    std::vector<int> result;
    if (step > 0) {
        for (int i = start; i < end; i += step)
            result.push_back(i);
    } else if (step < 0) {
        for (int i = start; i > end; i += step)
            result.push_back(i);
    }
    return result;
}

// ── identify_surfaces ─────────────────────────────────────────────────────────
// No Real needed — only accesses isol.stagIndex[] (always int[]).
template<typename IsolT, typename VsolT>
void identify_surfaces(const IsolT& isol, VsolT& vsol) {
    vsol.Is.clear();
    vsol.Is.push_back(range(isol.stagIndex[0], -1, -1));
    vsol.Is.push_back(range(isol.stagIndex[1], Ncoords));
    vsol.Is.push_back(range(Ncoords, Ncoords+Nwake));
}

// ── space_wake_nodes ─────────────────────────────────────────────────────────
// Real is in the parameter types so it is deduced at each call site.
template<typename Real, typename FoilT, typename WakeT>
void space_wake_nodes(const Real& wakeLength, const Real& firstPanelLength,
                      Real* wakeSpacing, const FoilT& foil, WakeT& wake) {

    int Nintervals = Nwake - 1;
    Real d = wakeLength / firstPanelLength;
    Real a = Nintervals*(Nintervals-1.0)*(Nintervals-2.0)/6.0;
    Real b = Nintervals*(Nintervals-1.0)/2.0;
    Real c = Nintervals - d;

    Real disc = std::max(b*b - 4.0*a*c, Real(0.0));
    Real r    = 1 + (-b + std::sqrt(disc)) / (2*a);

    // Newton-Raphson iterations
    Real R, R_r, dr;
    for (int k = 0; k < 10; ++k) {
        R   = std::pow(r, Nintervals) - 1 - d * (r - 1);
        R_r = Nintervals * std::pow(r, Nintervals - 1) - d;
        dr  = -R / R_r;
        if (std::abs(dr) < 1e-6) break;
        r -= R / R_r;
    }

    wakeSpacing[0] = 0.0;
    Real foilEndS = foil.s[Ncoords-1];
    wake.s[0] = foilEndS + 0.0;

    Real term = firstPanelLength;
    for (int i = 0; i < Nwake-1; ++i) {
        wakeSpacing[i+1] = wakeSpacing[i] + term;
        wake.s[i+1] = foilEndS + wakeSpacing[i+1];
        term *= r;
    }
}

// ── init_thermo ───────────────────────────────────────────────────────────────
// Real is deduced from oper.Vinf so no explicit Real template param is needed.
template<typename OperT, typename ParamT, typename GeomT>
void init_thermo(const OperT& oper, ParamT& param, const GeomT& geom) {

    using Real = std::decay_t<decltype(oper.Vinf)>;

    param.Vinf = oper.Vinf;
    param.muinf = oper.rho*oper.Vinf*geom.chord/oper.Re;
    if (oper.Ma > 0.0) {
        // Compressibility corrections
        Real gmi = param.gam - 1.0;
        param.KTb = std::sqrt(1.0 - oper.Ma*oper.Ma);
        param.KTl = (oper.Ma * oper.Ma) / std::pow(1.0 + param.KTb, 2);
        param.H0 = ((1.0 + 0.5 * gmi*oper.Ma*oper.Ma)*oper.Vinf*oper.Vinf) / (gmi*oper.Ma*oper.Ma);

        Real Tr = 1.0 - (0.5 * oper.Vinf*oper.Vinf) / param.H0;
        Real finf = std::pow(Tr, 1.5) * (1.0 + param.Tsrat) / (Tr + param.Tsrat);

        param.cps = (2.0 / (param.gam * oper.Ma*oper.Ma)) *
                    (std::pow((1.0 + 0.5*gmi*oper.Ma*oper.Ma) / (1.0 + 0.5 * gmi), param.gam / gmi) - 1.0);

        param.mu0 = param.muinf / finf;
        param.rho0 = oper.rho * std::pow(1.0 + 0.5 * gmi * oper.Ma * oper.Ma, 1.0 / gmi);
    } else {
        // Incompressible case
        param.mu0 = param.muinf;
        param.rho0 = oper.rho;
    }
}
