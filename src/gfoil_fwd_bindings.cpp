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
#include <stdexcept>

namespace py = pybind11;

static py::dict extract_obs(py::dict& d, const std::string& key) {
    return d;  // unused helper placeholder
}

py::dict run_forward_py(py::dict inp, py::object prev_jacobian = py::none()) {
    // ── geometry / aero inputs (runtime length nIn, floor NinMin) ────────────
    auto xlist = inp["xcoords"].cast<std::vector<double>>();
    auto ylist = inp["ycoords"].cast<std::vector<double>>();
    if (xlist.size() != ylist.size()) {
        throw std::invalid_argument(
            "xcoords and ycoords must be the same length; got " +
            std::to_string(xlist.size()) + " and " + std::to_string(ylist.size()));
    }
    const int nIn = static_cast<int>(xlist.size());
    if (nIn < NinMin) {
        throw std::invalid_argument(
            "input geometry needs at least " + std::to_string(NinMin) +
            " nodes for the cubic-spline re-panelling; got " + std::to_string(nIn));
    }
    std::vector<Real> inXcoords(nIn, Real(0.0));
    std::vector<Real> inYcoords(nIn, Real(0.0));
    for (int i = 0; i < nIn; ++i) {
        inXcoords[i] = xlist[i];
        inYcoords[i] = ylist[i];
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
        inXcoords.data(), inYcoords.data(), nIn,
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

// ── VALIDATION-ONLY entry point — not used by any production path ────────────
//
// Evaluates the Roger & Moreau radiation integral |I(omega, K2_bar)| with
// K2_bar as a FREE parameter, bypassing the Eq. 18 geometric gust selection
// K2_bar = k_bar*x2/S0 that TE_noise_outer_vec applies. That selection makes
// xi = beta*|x2|/S0 < 1 always, so production never reaches the subcritical
// branch except through the near-cutoff bridge; driving K2_bar independently is
// the only way to exercise the cut, as the paper's Figs. 9-11 do.
//
// The observer here is mid-span (S0 = sqrt(x^2 + beta^2 z^2)); K2_bar is
// supplied, not derived from it. Exposed as gfoil_cpp.amiet_kernel_I.
py::dict amiet_kernel_I_py(py::dict inp) {
    const double chord = inp["chord"].cast<double>();
    const double M     = inp["M"].cast<double>();
    const double Ue    = inp["Ue"].cast<double>();
    const double x     = inp["x"].cast<double>();
    const double z     = inp["z"].cast<double>();
    const bool bridged = inp["bridged"].cast<bool>();
    auto freqs  = inp["freqs_Hz"].cast<std::vector<double>>();
    auto K2_bar = inp["K2_bar"].cast<std::vector<double>>();

    if (freqs.size() != K2_bar.size())
        throw std::runtime_error("amiet_kernel_I: freqs_Hz and K2_bar must be "
                                 "the same length (they are paired per sample)");

    Real::getTape().reset();   // errFunc pushes statements; never evaluated here

    const double c0 = 340.0;
    const Real b     = Real(chord / 2.0);
    const Real U     = Real(M * c0);
    const Real beta  = std::sqrt(Real(1.0) - Real(M)*Real(M));
    const Real S0    = std::sqrt(Real(x)*Real(x) + beta*beta*Real(z)*Real(z));
    const Real alpha = U / (Real(0.7) * Real(Ue));

    const std::size_t N = freqs.size();
    std::vector<double> I_abs(N), xi_out(N), kappa_out(N), mu_out(N);
    std::vector<double> xi_a_out(N), xi_b_out(N), l_y_out(N);

    for (std::size_t i = 0; i < N; ++i) {
        Real omega   = Real(2.0 * M_PI * freqs[i]);
        Real K_bar   = (omega / U) * b;
        Real mu_bar  = K_bar * Real(M) / (beta*beta);
        Real K_1_bar = alpha * K_bar;
        Real C       = K_1_bar - mu_bar * (Real(x)/S0 - Real(M));
        Real xi      = std::abs(Real(K2_bar[i])) / (beta * mu_bar);

        Real Iabs = bridged
            ? Amiet_I_abs<Real>(xi, mu_bar, K_bar, K_1_bar, C, S0,
                                Real(x), Real(M), alpha)
            : Amiet_I_abs_raw<Real>(xi, mu_bar, K_bar, K_1_bar, C, S0,
                                    Real(x), Real(M), alpha);

        // Echo the bridge window so tests can check endpoint slope matching.
        Real kr = Real(AMIET_KAPPA_REG);
        if (mu_bar <= 2.0*AMIET_KAPPA_REG) kr = 0.5*mu_bar;
        Real t = kr / mu_bar;

        // Corcos length at the dimensional K₂ = K̄₂/b implied by this sample.
        l_y_out[i] = corcos_l_y<Real>(Real(K2_bar[i]) / b, omega, Real(1.47),
                                      Real(0.7) * Real(Ue)).getValue();

        I_abs[i]    = Iabs.getValue();
        xi_out[i]   = xi.getValue();
        mu_out[i]   = mu_bar.getValue();
        kappa_out[i] = (xi < 1.0)
            ? (mu_bar * std::sqrt(1.0 - xi*xi)).getValue()
            : -(mu_bar * std::sqrt(xi*xi - 1.0)).getValue();   // <0 flags kappa'
        xi_a_out[i] = std::sqrt(1.0 - t*t).getValue();
        xi_b_out[i] = std::sqrt(1.0 + t*t).getValue();
    }

    py::dict r;
    r["I_abs"]  = I_abs;
    r["xi"]     = xi_out;
    r["mu_bar"] = mu_out;
    r["kappa"]  = kappa_out;    // >0 supercritical kappa_bar, <0 subcritical -kappa'
    r["xi_a"]   = xi_a_out;
    r["xi_b"]   = xi_b_out;
    r["l_y"]    = l_y_out;      // R&M Eq. 19 at K_2 = K2_bar/b
    return r;
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
    m.def("amiet_kernel_I", &amiet_kernel_I_py,
          py::arg("input_dict"),
          "VALIDATION ONLY. |I(omega, K2_bar)| from the general Roger & Moreau "
          "radiation integral, with K2_bar a free parameter (bypasses the Eq. 18 "
          "geometric selection). Not used by any production path.");
}
