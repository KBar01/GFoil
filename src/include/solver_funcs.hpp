#pragma once

// Shared template implementations used by both GFoil_fwd_codi and
// GFoil_AD. Replaces duplicated non-template functions that previously
// existed in src/*.cpp and srcAD/include/main_func.hpp.
// Duck-typed on struct params (FoilT, IsolT, etc.) so it works with
// both the non-template fwd structs and the template AD structs.

#include <cmath>
#include <vector>

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
