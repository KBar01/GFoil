// Forward-pass pybind11 binding. Includes real_type.h (defines Real).
// Must NOT include real_type.hpp in the same TU.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "real_type.h"
#include "data_structs.h"
#include "restart_state.h"
#include "run_forward.h"
#include "gfoil_ad_impl.h"
#include "noise_run.hpp"
#include <string>
#include <vector>
#include <array>

namespace py = pybind11;

static py::dict extract_obs(py::dict& d, const std::string& key) {
    return d;  // unused helper placeholder
}

py::dict run_forward_py(py::dict inp, py::object prev_jacobian = py::none()) {
    // ── geometry / aero inputs ────────────────────────────────────────────────
    Real inXcoords[Nin] = {0};
    Real inYcoords[Nin] = {0};
    {
        auto xlist = inp["xcoords"].cast<std::vector<double>>();
        auto ylist = inp["ycoords"].cast<std::vector<double>>();
        for (int i = 0; i < Nin; ++i) {
            inXcoords[i] = xlist[i];
            inYcoords[i] = ylist[i];
        }
    }

    Real alphad      = inp["alpha_degrees"].cast<double>();
    Real Re          = inp["Re"].cast<double>();
    Real Ma          = inp["Ma"].cast<double>();
    Real rhoInf      = inp["rho"].cast<double>();
    Real kinViscInf  = inp["nu"].cast<double>();
    Real custChord   = inp["chord"].cast<double>();
    Real sampleTE    = inp["sampleTE"].cast<double>();
    Real sampleTE_hi = inp.contains("sampleTE_hi") ? inp["sampleTE_hi"].cast<double>()
                                                   : inp["sampleTE"].cast<double>();
    Real S           = inp["S"].cast<double>();
    Real Ncrit       = inp["ncrit"].cast<double>();
    Real Ufac        = inp["Ufac"].cast<double>();
    Real TEfac       = inp["TEfac"].cast<double>();
    int  doRestart   = inp["restart"].cast<int>();
    int  aWeighting  = inp.contains("aWeighting") ? inp["aWeighting"].cast<int>() : 0;
    bool verbose     = inp.contains("verbose")    ? inp["verbose"].cast<bool>() : false;
    double f_min     = inp.contains("f_min")      ? inp["f_min"].cast<double>() : 200.0;
    double f_max     = inp.contains("f_max")      ? inp["f_max"].cast<double>() : 20000.0;
    double xft_lower = inp.contains("bottrans")   ? inp["bottrans"].cast<double>() : 1.0;
    double xft_upper = inp.contains("toptrans")   ? inp["toptrans"].cast<double>() : 1.0;
    double rtol      = inp.contains("rtol")       ? inp["rtol"].cast<double>() : 1e-6;
    std::string model = inp["model"].cast<std::string>();

    // ── observer arrays ───────────────────────────────────────────────────────
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
    std::vector<Real> obsX(nObs), obsY(nObs), obsZ(nObs);
    for (int i = 0; i < nObs; ++i) {
        obsX[i] = obsX_d[i];
        obsY[i] = obsY_d[i];
        obsZ[i] = obsZ_d[i];
    }

    // ── warm-start from previous result (pybind11 continuation path) ─────────
    RestartState warmStartState;
    const RestartState* warmStartPtr = nullptr;
    if (!prev_jacobian.is_none()) {
        py::dict jac = prev_jacobian.cast<py::dict>();
        auto states_py = jac["states"].cast<std::vector<double>>();
        auto turb_py   = jac["turb"].cast<std::vector<int>>();
        warmStartState.states.assign(states_py.begin(), states_py.end());
        warmStartState.turb.assign(turb_py.begin(), turb_py.end());
        // Donor's converged stagnation bracket — used to seed stagpoint_move on
        // warm entry so it reproduces the donor configuration instead of landing
        // one node off the inviscid seed. Optional for backward compatibility.
        if (jac.contains("stag")) {
            auto stag_py = jac["stag"].cast<std::vector<int>>();
            warmStartState.stag.assign(stag_py.begin(), stag_py.end());
        }
        warmStartPtr = &warmStartState;
    }

    // ── run solver ────────────────────────────────────────────────────────────
    RestartState rst;
    ForwardResult fwd;
    bool converged = runCode(
        static_cast<bool>(doRestart),
        Ncrit, Ufac, TEfac, custChord,
        inXcoords, inYcoords,
        alphad, Re, Ma, rhoInf, kinViscInf,
        model, sampleTE, sampleTE_hi,
        obsX.data(), obsY.data(), obsZ.data(), nObs,
        S,
        &rst, &fwd,
        warmStartPtr,
        aWeighting,
        verbose,
        f_min,
        f_max,
        xft_lower,
        xft_upper,
        rtol);

    // ── pack result ───────────────────────────────────────────────────────────
    py::dict result;
    result["conv"]              = converged ? 1 : 0;
    result["failure_mode"]      = fwd.failure_mode;
    result["newton_iterations"] = fwd.newton_iterations;
    if (converged) {
        result["CL"]    = fwd.CL;
        result["CD"]    = fwd.CD;
        result["CM"]    = fwd.CM;
        result["OASPL"] = fwd.OASPL;

        py::dict jac;
        jac["states"] = rst.states;
        jac["turb"]   = rst.turb;
        jac["stag"]   = rst.stag;
        jac["RVvals"] = rst.RVvals;
        jac["RVrows"] = rst.RVrows;
        jac["RVcols"] = rst.RVcols;
        jac["RVnz"]   = rst.RVnz;
        result["jacobian"] = jac;

        // ── verbose per-node and acoustic data ────────────────────────────────
        if (!fwd.innerFoilX.empty()) {
            result["innerFoilX"]  = fwd.innerFoilX;
            result["innerFoilY"]  = fwd.innerFoilY;
            result["Cp_dist"]     = fwd.Cp;
            result["delta_star"]  = fwd.delta_star;
            result["theta"]       = fwd.theta;
            result["tau_wall"]    = fwd.tau_wall;
            result["tau_max"]     = fwd.tau_max;
            result["Ue"]          = fwd.Ue;
            result["dpdx"]        = fwd.dpdx;
            result["is_turb"]     = fwd.is_turb;
            result["topTransX"]   = fwd.topTransX;
            result["botTransX"]   = fwd.botTransX;
            result["BL_top"]      = fwd.BL_top;
            result["BL_bot"]      = fwd.BL_bot;
            result["freq_Hz"]     = fwd.freq_Hz;
            result["WPS_upper"]   = fwd.WPS_upper;
            result["WPS_lower"]   = fwd.WPS_lower;
            result["FF_spectra"]  = fwd.FF_spectra;
            result["nObs"]        = fwd.nObs;
            result["OASPL_perObs"]   = fwd.OASPL_perObs;
            result["obsXYZ_TElocal"] = fwd.obsXYZ_TElocal;
        }
    }
    return result;
}

// ── acoustics-only entry point (no aero solve, never AD'd) ───────────────────
py::dict run_noise_py(py::dict inp) {
    double alphaDeg = inp["alphaDeg"].cast<double>();
    double Re       = inp["Re"].cast<double>();
    double rho      = inp["rho"].cast<double>();
    double nu       = inp["nu"].cast<double>();
    double Ma       = inp["Ma"].cast<double>();
    double chord    = inp["chord"].cast<double>();
    double span     = inp["span"].cast<double>();
    std::string model = inp["model"].cast<std::string>();

    // ── observer arrays (accept list-or-scalar, like run_forward_py) ──────────
    std::vector<double> obsX, obsY, obsZ;
    auto xval = inp["X"];
    if (py::isinstance<py::list>(xval) || py::isinstance<py::sequence>(xval)) {
        obsX = xval.cast<std::vector<double>>();
        obsY = inp["Y"].cast<std::vector<double>>();
        obsZ = inp["Z"].cast<std::vector<double>>();
    } else {
        obsX = { xval.cast<double>() };
        obsY = { inp["Y"].cast<double>() };
        obsZ = { inp["Z"].cast<double>() };
    }

    // ── BL state pairs [upper, lower] ─────────────────────────────────────────
    auto theta     = inp["theta"].cast<std::array<double,2>>();
    auto deltaStar = inp["deltaStar"].cast<std::array<double,2>>();
    auto tauMax    = inp["tauMax"].cast<std::array<double,2>>();
    auto Ue        = inp["Ue"].cast<std::array<double,2>>();
    auto dpdx      = inp["dpdx"].cast<std::array<double,2>>();
    auto tauWall   = inp["tauWall"].cast<std::array<double,2>>();
    auto delta99   = inp["delta99"].cast<std::array<double,2>>();

    auto freqs_Hz  = inp["freqs_Hz"].cast<std::vector<double>>();

    // ── optional custom WPS, shape (N,2) columns [upper, lower] ───────────────
    bool has_custom = inp.contains("custom_WPS");
    std::vector<double> custom_upper, custom_lower;
    if (has_custom) {
        auto cw = inp["custom_WPS"].cast<std::vector<std::array<double,2>>>();
        custom_upper.resize(cw.size());
        custom_lower.resize(cw.size());
        for (std::size_t i = 0; i < cw.size(); ++i) {
            custom_upper[i] = cw[i][0];
            custom_lower[i] = cw[i][1];
        }
    }

    NoiseRunResult r = noise_run_cpp<Real>(
        alphaDeg, Re, rho, nu, Ma, chord,
        obsX, obsY, obsZ, span,
        theta, deltaStar, tauMax, Ue, dpdx, tauWall, delta99,
        freqs_Hz, model,
        has_custom, custom_upper, custom_lower);

    py::dict result;
    result["freqs_Hz"]       = r.freqs_Hz;
    result["WPS_upper"]      = r.WPS_upper;
    result["WPS_lower"]      = r.WPS_lower;
    result["FF_spectra"]     = r.FF_spectra;        // flat nObs*N, row-major
    result["obsXYZ_TElocal"] = r.obsXYZ_TElocal;    // flat nObs*3
    result["nObs"]           = static_cast<int>(obsX.size());
    result["N"]              = static_cast<int>(freqs_Hz.size());
    return result;
}

PYBIND11_MODULE(gfoil_cpp, m) {
    m.def("run_forward", &run_forward_py,
          py::arg("input_dict"),
          py::arg("prev_jacobian") = py::none(),
          "Run forward aero+acoustic solver. Returns dict with CL/CD/CM/OASPL and jacobian.");
    m.def("run_AD", &run_AD_py,
          "Run AD solver given input dict and jacobian from run_forward.");
    m.def("noise_run", &run_noise_py,
          py::arg("input_dict"),
          "Acoustics-only forward run (no aero solve, never AD'd). Returns raw "
          "linear Pa^2/omega WPS and far-field spectra for arbitrary-length freqs.");
}
