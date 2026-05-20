// AD-pass pybind11 binding. Includes real_type.hpp — must NOT include real_type.h
// (both define norm2; including both in one TU causes a redefinition error).

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// real_type.hpp FIRST — defines macros (Nin, Ncoords, RVdimension, Nsound …)
// and norm2 as a template. DO NOT include real_type.h after this point.
#include "real_type.hpp"   // from srcAD/include/ via include path

#include "codi.hpp"
#include "ADfuncs.hpp"     // from srcAD/include/

#include "gfoil_ad_impl.h"

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

namespace py = pybind11;

using MatrixSparse = Eigen::SparseMatrix<double>;
using VectorD      = Eigen::Matrix<double, Eigen::Dynamic, 1>;

// solve_sys is defined in srcAD/main.cpp which is NOT compiled into the module.
// Replicate it here (identical logic).
static void solve_sys_ad(int* R_V_cols, int* R_V_rows, double* RV_vals, const int RVlatest,
    double (&dCLdStates)[RVdimension],
    double (&dCDdStates)[RVdimension],
    double (&dOASPLdStates)[RVdimension],
    double (&adlambdaCL)[RVdimension],
    double (&adlambdaCD)[RVdimension],
    double (&adlambdaOASPL)[RVdimension])
{
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(RVlatest);
    for (int k = 0; k < RVlatest; ++k)
        triplets.emplace_back(R_V_rows[k], R_V_cols[k], RV_vals[k]);

    MatrixSparse A(RVdimension, RVdimension);
    A.setFromTriplets(triplets.begin(), triplets.end());

    VectorD rhs1(RVdimension), rhs2(RVdimension), rhs3(RVdimension);
    for (int i = 0; i < RVdimension; ++i) {
        rhs1(i) = dCLdStates[i];
        rhs2(i) = dCDdStates[i];
        rhs3(i) = dOASPLdStates[i];
    }

    Eigen::SparseLU<MatrixSparse> solver;
    solver.compute(A);
    VectorD sol1 = solver.solve(rhs1);
    VectorD sol2 = solver.solve(rhs2);
    VectorD sol3 = solver.solve(rhs3);

    for (int i = 0; i < RVdimension; ++i) {
        adlambdaCL[i]    = sol1(i);
        adlambdaCD[i]    = sol2(i);
        adlambdaOASPL[i] = sol3(i);
    }
}

py::dict run_AD_py(py::dict inp, py::dict jacobian) {
    using RealVec2 = codi::RealReverseVec<2>;
    using RealRev  = codi::RealReverseVec<3>;
    using Realfwd  = codi::RealForward;

    // ── unpack jacobian ───────────────────────────────────────────────────────
    auto states_py = jacobian["states"].cast<std::vector<double>>();
    auto turb_py   = jacobian["turb"].cast<std::vector<int>>();
    auto stag_py   = jacobian["stag"].cast<std::vector<int>>();
    auto RVvals_py = jacobian["RVvals"].cast<std::vector<double>>();
    auto RVrows_py = jacobian["RVrows"].cast<std::vector<int>>();
    auto RVcols_py = jacobian["RVcols"].cast<std::vector<int>>();
    int  RVnz      = jacobian["RVnz"].cast<int>();

    // Static arrays to avoid stack overflow (~900 KB each)
    static RealVec2 states[RVdimension];
    static double   states_d[RVdimension];
    for (int i = 0; i < RVdimension; ++i) {
        states[i]   = states_py[i];
        states_d[i] = states_py[i];
    }

    int turb[Ncoords + Nwake];
    for (int i = 0; i < Ncoords + Nwake; ++i)
        turb[i] = turb_py[i];

    int currStag[2] = {stag_py[0], stag_py[1]};

    // Rows↔cols swap preserved from srcAD/main.cpp read path
    static double dRdU_vals[119700];
    static int    dRdU_rows[119700];
    static int    dRdU_cols[119700];
    for (int i = 0; i < RVnz; ++i) {
        dRdU_vals[i] = RVvals_py[i];
        dRdU_rows[i] = RVcols_py[i];   // intentional swap (mirrors restart.json read)
        dRdU_cols[i] = RVrows_py[i];
    }

    // ── unpack input ─────────────────────────────────────────────────────────
    double inXcoords_d[Nin] = {0};
    RealVec2 inYcoords_2[Nin] = {};
    RealRev  inYcoords_Rev[Nin] = {};
    {
        auto xlist = inp["xcoords"].cast<std::vector<double>>();
        auto ylist = inp["ycoords"].cast<std::vector<double>>();
        for (int i = 0; i < Nin; ++i) {
            inXcoords_d[i]    = xlist[i];
            inYcoords_2[i]    = ylist[i];
            inYcoords_Rev[i]  = ylist[i];
        }
    }

    RealVec2 targetAlphaDeg = inp["alpha_degrees"].cast<double>();
    RealVec2 Re             = inp["Re"].cast<double>();
    RealVec2 Ma             = inp["Ma"].cast<double>();
    RealVec2 rhoInf         = inp["rho"].cast<double>();
    RealVec2 nuInf          = inp["nu"].cast<double>();
    RealVec2 custChord      = inp["chord"].cast<double>();
    RealVec2 sampleTE       = inp["sampleTE"].cast<double>();
    const RealVec2 S        = inp["S"].cast<double>();
    const RealVec2 Ncrit    = inp["ncrit"].cast<double>();
    const RealVec2 Ufac     = inp["Ufac"].cast<double>();
    const RealVec2 TEfac    = inp["TEfac"].cast<double>();
    std::string model       = inp["model"].cast<std::string>();
    int aWeighting          = inp.contains("aWeighting") ? inp["aWeighting"].cast<int>() : 0;

    RealRev targetAlphaDeg_r = inp["alpha_degrees"].cast<double>();
    RealRev Re_r             = inp["Re"].cast<double>();
    RealRev Ma_r             = inp["Ma"].cast<double>();
    RealRev rhoInf_r         = inp["rho"].cast<double>();
    const RealRev Ncrit_r    = inp["ncrit"].cast<double>();
    const RealRev Ufac_r     = inp["Ufac"].cast<double>();
    const RealRev TEfac_r    = inp["TEfac"].cast<double>();

    // observer arrays — both Vec2 and Rev types
    std::vector<double> obsX_d, obsY_d, obsZ_d;
    auto xval = inp["X"];
    if (py::isinstance<py::list>(xval) || py::isinstance<py::sequence>(xval)) {
        obsX_d = xval.cast<std::vector<double>>();
        obsY_d = inp["Y"].cast<std::vector<double>>();
        obsZ_d = inp["Z"].cast<std::vector<double>>();
    } else {
        obsX_d = { xval.cast<double>() };
        obsY_d = { inp["Y"].cast<double>() };
        obsZ_d = { inp["Z"].cast<double>() };
    }
    int nObs = static_cast<int>(obsX_d.size());
    std::vector<RealVec2> obsX_2(nObs), obsY_2(nObs), obsZ_2(nObs);
    for (int i = 0; i < nObs; ++i) {
        obsX_2[i] = obsX_d[i];
        obsY_2[i] = obsY_d[i];
        obsZ_2[i] = obsZ_d[i];
    }

    // ── kinematic chain term (CD via Realfwd, same as srcAD/main.cpp) ────────
    const double theta_d = states_py[RVdimension - 4];
    const double disp_d  = states_py[RVdimension - 3];
    const double ue_d    = states_py[RVdimension - 1];

    Realfwd exponent = 2.5 + disp_d / (2.0 * theta_d);
    Realfwd ue_pow   = std::pow(ue_d, exponent);
    Realfwd dCdcmom  = ue_pow * (2.0 * theta_d - disp_d * std::log(ue_d)) / theta_d;
    Realfwd Cdddisp  = std::log(ue_d) * ue_pow;
    Realfwd Cddue    = (2.0 * theta_d) * exponent * std::pow(ue_d, exponent - 1.0);

    // ── AD pass 1: partialOutputspartialInputs ────────────────────────────────
    double d_CL_d_y[Nin]    = {0};
    double d_OASPL_d_y[Nin] = {0};
    double d_CL_dalpha = 0.0, d_OASPL_dalpha = 0.0;
    static double d_CL_d_States[RVdimension]    = {0};
    static double d_CD_d_States[RVdimension]    = {0};
    static double d_OASPL_d_States[RVdimension] = {0};

    partialOutputspartialInputs<RealVec2>(
        Ncrit, Ufac, TEfac, custChord, inXcoords_d, Re, Ma, rhoInf, nuInf,
        model, sampleTE,
        obsX_2.data(), obsY_2.data(), obsZ_2.data(), nObs,
        S, inYcoords_2, targetAlphaDeg, states, turb,
        d_CL_d_y, d_OASPL_d_y, d_CL_dalpha, d_OASPL_dalpha,
        d_CL_d_States, d_OASPL_d_States,
        aWeighting);

    // ── adjoint solve ─────────────────────────────────────────────────────────
    d_CD_d_States[RVdimension - 4] = dCdcmom.getValue();
    d_CD_d_States[RVdimension - 3] = Cdddisp.getValue();
    d_CD_d_States[RVdimension - 1] = Cddue.getValue();

    static double adlambda_CL[RVdimension]    = {0.0};
    static double adlambda_CD[RVdimension]    = {0.0};
    static double adlambda_OASPL[RVdimension] = {0.0};

    solve_sys_ad(dRdU_cols, dRdU_rows, dRdU_vals, RVnz,
        d_CL_d_States, d_CD_d_States, d_OASPL_d_States,
        adlambda_CL, adlambda_CD, adlambda_OASPL);

    // ── AD pass 2: partialRpartialx ───────────────────────────────────────────
    static double dgdy_CL[Nin]    = {0};
    static double dgdy_CD[Nin]    = {0};
    static double dgdy_OASPL[Nin] = {0};
    double dgdalpha_CL = 0, dgdalpha_CD = 0, dgdalpha_OASPL = 0;

    partialRpartialx<RealRev>(
        Ncrit_r, Ufac_r, TEfac_r, inXcoords_d, Re_r, Ma_r, rhoInf_r, currStag,
        adlambda_CL, adlambda_CD, adlambda_OASPL,
        inYcoords_Rev, targetAlphaDeg_r, states_d, turb,
        dgdy_CL, dgdy_CD, dgdy_OASPL, dgdalpha_CL, dgdalpha_CD, dgdalpha_OASPL);

    // ── total derivatives ─────────────────────────────────────────────────────
    std::vector<double> totalCL(Nin), totalCD(Nin), totalOASPL(Nin);
    for (int i = 0; i < Nin; ++i) {
        totalCL[i]    = d_CL_d_y[i]    + dgdy_CL[i];
        totalCD[i]    = dgdy_CD[i];
        totalOASPL[i] = d_OASPL_d_y[i] + dgdy_OASPL[i];
    }

    py::dict grads;
    grads["dCL_dy"]        = totalCL;
    grads["dCL_dalpha"]    = d_CL_dalpha + dgdalpha_CL;
    grads["dCD_dy"]        = totalCD;
    grads["dCD_dalpha"]    = dgdalpha_CD;
    grads["dOASPL_dy"]     = totalOASPL;
    grads["dOASPL_dalpha"] = d_OASPL_dalpha + dgdalpha_OASPL;
    return grads;
}
