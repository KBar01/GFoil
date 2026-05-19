#pragma once

// Shared template implementations used by both GFoil_fwd_codi and
// GFoil_AD. Replaces duplicated non-template functions that previously
// existed in src/*.cpp and srcAD/include/main_func.hpp.
// Duck-typed on struct params (FoilT, IsolT, etc.) so it works with
// both the non-template fwd structs and the template AD structs.
//
// Functions defined here:
//   build_wake_impl<FoilT,GeomT,OperT,IsolcT,WakeT>
//   calc_force<OperT,GeomT,ParamT,FoilT,GlobT,PostT>
//   stagpoint_find_impl<bool,GammaT,VarT,FoilT,WakeT>
//   rebuild_ue_m<FoilT,WakeT,IsolT,VsolT>
//   set_wake_gap<FoilT,IsolT,VsolT>
//   range(int,int,int)
//   identify_surfaces<IsolT,VsolT>
//   space_wake_nodes<Real,FoilT,WakeT>
//   init_thermo<OperT,ParamT,GeomT>

#include <cmath>
#include <vector>

// ── build_wake_impl ──────────────────────────────────────────────────────────
// Duck-typed on IsolcT — both Isol (fwd) and Isolc<Real> (AD) provide
// .gammas[] and .uewi[]. Real deduced from op.Vinf.
template<typename FoilT, typename GeomT, typename OperT,
         typename IsolcT, typename WakeT>
void build_wake_impl(const FoilT& foil, const GeomT& geom,
                     const OperT& op, IsolcT& isol, WakeT& wake) {

    using Real = std::decay_t<decltype(op.Vinf)>;

    Real firstPanelSize = 0.5*(foil.s[1]-foil.s[0] + foil.s[Ncoords-1]-foil.s[Ncoords-2]);
    Real wakeLength = geom.wakelen*geom.chord;
    Real wakePanelSizes[Nwake];
    space_wake_nodes(wakeLength,firstPanelSize,wakePanelSizes,foil,wake);

    Real midpointTE[2];
    midpointTE[0] = 0.5*(foil.x[colMajorIndex(0,0,2)] + foil.x[colMajorIndex(0,Ncoords-1,2)]);
    midpointTE[1] = 0.5*(foil.x[colMajorIndex(1,0,2)] + foil.x[colMajorIndex(1,Ncoords-1,2)]);

    Real normTE[2];
    Real tangTE[2];
    normTE[0] = foil.x[colMajorIndex(0,Ncoords-1,2)] - foil.x[colMajorIndex(0,0,2)];
    normTE[1] = foil.x[colMajorIndex(1,Ncoords-1,2)] - foil.x[colMajorIndex(1,0,2)];

    tangTE[0] = normTE[1];
    tangTE[1] = -normTE[0];

    wake.x[0] = midpointTE[0] + (1.0e-5)*tangTE[0]*geom.chord;
    wake.x[1] = midpointTE[1] + (1.0e-5)*tangTE[1]*geom.chord;

    Real v1[2],v2[2],v1Store[2];
    Real norm,wakeLengthsDiff;
    for (int i=0; i<Nwake-1; ++i) {

        v1Store[0]=0.0; v1Store[1]=0.0;
        v2[0]=0.0;      v2[1]=0.0;
        inviscid_velocity(foil,isol.gammas,op.Vinf,op.alpha,
                          wake.x[colMajorIndex(0,i,2)],wake.x[colMajorIndex(1,i,2)],v1Store);

        norm = norm2(v1Store);
        v1[0] = v1Store[0]/norm;
        v1[1] = v1Store[1]/norm;
        wake.t[colMajorIndex(0,i,2)] = v1[0];
        wake.t[colMajorIndex(1,i,2)] = v1[1];

        isol.uewi[i] = (v1Store[0] * v1[0]) + (v1Store[1] * v1[1]);

        wakeLengthsDiff = wakePanelSizes[i+1]-wakePanelSizes[i];
        wake.x[colMajorIndex(0,i+1,2)] = wake.x[colMajorIndex(0,i,2)] + wakeLengthsDiff*v1[0];
        wake.x[colMajorIndex(1,i+1,2)] = wake.x[colMajorIndex(1,i,2)] + wakeLengthsDiff*v1[1];

        inviscid_velocity(foil,isol.gammas,op.Vinf,op.alpha,
                          wake.x[colMajorIndex(0,i+1,2)],wake.x[colMajorIndex(1,i+1,2)],v2);
        norm = norm2(v2);
        v2[0] /= norm; v2[1] /= norm;
        wake.t[colMajorIndex(0,i+1,2)] = v2[0];
        wake.t[colMajorIndex(1,i+1,2)] = v2[1];

        wake.x[colMajorIndex(0,i+1,2)] = wake.x[colMajorIndex(0,i,2)] + wakeLengthsDiff*0.5*(v1[0]+v2[0]);
        wake.x[colMajorIndex(1,i+1,2)] = wake.x[colMajorIndex(1,i,2)] + wakeLengthsDiff*0.5*(v1[1]+v2[1]);
    }

    v1Store[0]=0.0; v1Store[1]=0.0;
    inviscid_velocity(foil,isol.gammas,op.Vinf,op.alpha,
                      wake.x[colMajorIndex(0,Nwake-1,2)],wake.x[colMajorIndex(1,Nwake-1,2)],v1Store);
    norm = norm2(v1Store);
    v1[0] = v1Store[0]/norm;
    v1[1] = v1Store[1]/norm;
    wake.t[colMajorIndex(0,Nwake-1,2)] = v1[0];
    wake.t[colMajorIndex(1,Nwake-1,2)] = v1[1];

    isol.uewi[Nwake-1] = (v1Store[0] * v1[0]) + (v1Store[1] * v1[1]);
}

// ── calc_force ───────────────────────────────────────────────────────────────
// Duck-typed on all struct params. Real deduced from glob.U element type.
// The fwd signature has an extra unused isol param (dropped here).
template<typename OperT, typename GeomT, typename ParamT,
         typename FoilT, typename GlobT, typename PostT>
void calc_force(const OperT& op, const GeomT& geom, const ParamT& par,
                const FoilT& foil, const GlobT& glob, PostT& post) {

    using Real = std::decay_t<decltype(glob.U[0])>;

    Real qinf = 0.5*op.rho *(op.Vinf * op.Vinf);
    constexpr int Nsys = Ncoords+Nwake;

    Real ue[Nsys]={0};
    for (int i=0;i<Nsys;++i){ue[i] = glob.U[colMajorIndex(3,i,4)];}
    get_cp(post,op,par,ue);

    Real cos_alpha = std::cos(op.alpha);
    Real sin_alpha = std::sin(op.alpha);

    for (int i0 = 1; i0 <= Ncoords; ++i0) {

        int i  = (i0 == Ncoords) ? 0        : i0;
        int ip = (i0 == Ncoords) ? Ncoords-1 : i0-1;

        const Real x1[2] = {foil.x[colMajorIndex(0,ip,2)],foil.x[colMajorIndex(1,ip,2)]};
        const Real x2[2] = {foil.x[colMajorIndex(0,i ,2)],foil.x[colMajorIndex(1,i ,2)]};

        Real dxv[2] = {x2[0] - x1[0], x2[1] - x1[1]};
        Real dx1[2] = {x1[0] - geom.xref[0], x1[1] - geom.xref[1]};
        Real dx2[2] = {x2[0] - geom.xref[0], x2[1] - geom.xref[1]};

        Real dx1nds = dxv[0] * dx1[0] + dxv[1] * dx1[1];
        Real dx2nds = dxv[0] * dx2[0] + dxv[1] * dx2[1];

        Real dx = -dxv[0] * cos_alpha - dxv[1] * sin_alpha;
        Real dz =  dxv[1] * cos_alpha - dxv[0] * sin_alpha;

        Real cp1 = post.cp[ip], cp2 = post.cp[i];
        Real cpbar = 0.5 * (cp1 + cp2);

        post.cl += dx * cpbar;
        post.cm += cp1 * dx1nds / 3.0 + cp1 * dx2nds / 6.0
                 + cp2 * dx1nds / 6.0 + cp2 * dx2nds / 3.0;
    }

    post.cl /= geom.chord;
    post.cm /= (geom.chord*geom.chord);

    int iw = Ncoords+Nwake-1;
    const Real* U = &glob.U[colMajorIndex(0,iw,4)];

    Real H, H_U[4];
    H = get_H(U[0],U[1],H_U);

    Real uk, uk_ue;
    uk = get_uk(U[3],par,uk_ue);

    post.cd = 2.0 * U[0] * pow(uk / op.Vinf, (5.0 + H) / 2.0);
}

// ── stagpoint_find_impl ───────────────────────────────────────────────────────
// GammaT provides .gammas[] (Isol in fwd, Isolc<Real> in AD).
// VarT    receives .stagIndex[], .stagArcLocation, etc. (Isol in fwd, Isolv<Real> in AD).
// compute_sstag_g is a compile-time bool: true for fwd (Isol has sstag_g[]),
// false for AD (Isolv does not). if constexpr prevents the field access from
// being instantiated for the AD path.
template<bool compute_sstag_g, typename GammaT, typename VarT,
         typename FoilT, typename WakeT>
void stagpoint_find_impl(const GammaT& gammaStruct, VarT& varStruct,
                         const FoilT& foil, const WakeT& wake) {
    int j = 0;
    for (; j < Ncoords; ++j) {
        if (gammaStruct.gammas[j] > 0) break;
    }
    int I[2] = {j-1, j};
    varStruct.stagIndex[0] = I[0];
    varStruct.stagIndex[1] = I[1];

    auto G0 = gammaStruct.gammas[I[0]];
    auto G1 = gammaStruct.gammas[I[1]];
    auto S0 = foil.s[I[0]];
    auto S1 = foil.s[I[1]];

    auto den = G1 - G0;
    auto w1  = G1 / den;
    auto w2  = -G0 / den;

    varStruct.stagArcLocation = w1 * S0 + w2 * S1;

    varStruct.stagXLocation[0] = foil.x[colMajorIndex(0,j-1,2)] * w1 + foil.x[colMajorIndex(1,j,2)] * w2;
    varStruct.stagXLocation[1] = foil.x[colMajorIndex(1,j-1,2)] * w1 + foil.x[colMajorIndex(1,j,2)] * w2;

    auto st_g1 = G1 * (S0 - S1) / (den * den);
    if constexpr (compute_sstag_g) {
        // sstag_g[] only exists in the fwd Isol struct, not in AD's Isolv
        varStruct.sstag_g[0] =  st_g1;
        varStruct.sstag_g[1] = -st_g1;
    }

    for (int i = j; i < Ncoords; ++i)
        varStruct.edgeVelSign[i] = 1.0;

    for (int i = 0; i < Ncoords; ++i)
        varStruct.distFromStag[i] = std::abs(foil.s[i] - varStruct.stagArcLocation);

    for (int i = 0; i < Nwake; ++i)
        varStruct.distFromStag[Ncoords + i] = wake.s[i] - varStruct.stagArcLocation;
}

// ── rebuild_ue_m ─────────────────────────────────────────────────────────────
// Real deduced from isol.edgeVelSign[0] element type.
template<typename FoilT, typename WakeT, typename IsolT, typename VsolT>
void rebuild_ue_m(const FoilT& foil, const WakeT& wake,
                  const IsolT& isol, VsolT& vsol, bool realloc) {

    using Real = std::decay_t<decltype(foil.s[0])>;  // edgeVelSign is int[], deduce from foil.s

    Real sigma_m[2*(Ncoords-1)];

    if (realloc) {
        for (int i = 0; i < (Ncoords+Nwake)*(Ncoords+Nwake); ++i)
            vsol.ue_m[i] = 0;
    }

    for (int i = 0; i < Ncoords-1; ++i) {
        Real ds = foil.s[i+1] - foil.s[i];
        int ind = 2*i;
        sigma_m[ind]   = -1.0*isol.edgeVelSign[i]   / ds;
        sigma_m[ind+1] =  1.0*isol.edgeVelSign[i+1] / ds;
    }

    for (int row = 0; row < Ncoords+Nwake; ++row)
        vsol.ue_m[row] = vsol.ue_sigma[row] * sigma_m[0];

    int ue_mInt = Ncoords+Nwake;
    int sigma_mIndex = 1;
    for (int col = 1; col <= Ncoords-2; ++col) {
        Real sigma_mValue = sigma_m[sigma_mIndex];
        for (int row = 0; row < Ncoords+Nwake; ++row)
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,col-1,Ncoords+Nwake)] * sigma_mValue;
        sigma_mValue = sigma_m[sigma_mIndex+1];
        for (int row = 0; row < Ncoords+Nwake; ++row)
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,col,Ncoords+Nwake)] * sigma_mValue;
        sigma_mIndex += 2;
        ue_mInt += (Ncoords+Nwake);
    }

    for (int row = 0; row < Ncoords+Nwake; ++row)
        vsol.ue_m[colMajorIndex(row,Ncoords-1,Ncoords+Nwake)] =
            vsol.ue_sigma[colMajorIndex(row,Ncoords-2,Ncoords+Nwake)] * sigma_m[2*(Ncoords-1)-1];

    Real sigma_mWake[2*(Nwake-1)];
    for (int i = 0; i < Nwake-1; ++i) {
        int ind = 2*i;
        Real ds = wake.s[i+1] - wake.s[i];
        sigma_mWake[ind]   = -1.0 / ds;
        sigma_mWake[ind+1] =  1.0 / ds;
    }

    for (int row = 0; row < Ncoords+Nwake; ++row)
        vsol.ue_m[colMajorIndex(row,Ncoords,Ncoords+Nwake)] =
            vsol.ue_sigma[colMajorIndex(row,Ncoords-1,Ncoords+Nwake)] * sigma_mWake[0];

    ue_mInt = (Ncoords+Nwake)*(Ncoords) + Ncoords+Nwake;
    sigma_mIndex = 1;
    for (int col = 1; col <= Nwake-2; ++col) {
        Real sigma_mValue = sigma_mWake[sigma_mIndex];
        for (int row = 0; row < Ncoords+Nwake; ++row)
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,Ncoords+col-2,Ncoords+Nwake)] * sigma_mValue;
        sigma_mValue = sigma_mWake[sigma_mIndex+1];
        for (int row = 0; row < Ncoords+Nwake; ++row)
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,Ncoords+col-1,Ncoords+Nwake)] * sigma_mValue;
        sigma_mIndex += 2;
        ue_mInt += (Ncoords+Nwake);
    }

    for (int row = 0; row < Ncoords+Nwake; ++row)
        vsol.ue_m[(Ncoords+Nwake)*(Ncoords+Nwake-1)+row] =
            vsol.ue_sigma[colMajorIndex(row,Ncoords+Nwake-3,Ncoords+Nwake)] * sigma_mWake[2*(Nwake-1)-1];

    for (int row = 0; row < Ncoords; ++row) {
        if (isol.edgeVelSign[row] == 1.0) continue;
        for (int col = 0; col < Ncoords+Nwake; ++col)
            vsol.ue_m[colMajorIndex(row,col,Ncoords+Nwake)] *= -1.0;
    }
}

// ── set_wake_gap ─────────────────────────────────────────────────────────────
// Accesses foil.te.hTE, foil.te.dtdx, isol.distFromStag[], vsol.wgap[].
template<typename FoilT, typename IsolT, typename VsolT>
void set_wake_gap(const FoilT& foil, const IsolT& isol, VsolT& vsol) {

    using Real = std::decay_t<decltype(foil.te.hTE)>;

    Real lengthScaleFactr = 2.5;
    Real dtdx = std::min(std::max(foil.te.dtdx, -3.0/lengthScaleFactr), 3.0/lengthScaleFactr);
    Real Lw = lengthScaleFactr * foil.te.hTE;

    Real xib;
    for (int i = 0; i < Nwake; ++i) {
        xib = (isol.distFromStag[Ncoords+i] - isol.distFromStag[Ncoords]) / Lw;
        if (xib <= 1.0) {
            vsol.wgap[i] = foil.te.hTE * (1 + (2 + lengthScaleFactr*dtdx)*xib) * (1-xib) * (1-xib);
        }
    }
}

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
