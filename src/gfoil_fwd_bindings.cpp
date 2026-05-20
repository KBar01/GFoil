// Forward-pass pybind11 binding. Includes real_type.h (defines Real).
// Must NOT include real_type.hpp in the same TU.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "real_type.h"
#include "data_structs.h"
#include "restart_state.h"
#include "run_forward.h"
#include "gfoil_ad_impl.h"
#include <string>
#include <vector>

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
    Real S           = inp["S"].cast<double>();
    Real Ncrit       = inp["ncrit"].cast<double>();
    Real Ufac        = inp["Ufac"].cast<double>();
    Real TEfac       = inp["TEfac"].cast<double>();
    int  doRestart   = inp["restart"].cast<int>();
    int  doCps       = inp["returnData"].cast<int>();
    int  aWeighting  = inp.contains("aWeighting") ? inp["aWeighting"].cast<int>() : 0;
    Real ncrithyst   = inp.contains("ncrithyst")  ? Real(inp["ncrithyst"].cast<double>()) : Real(0.2);
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
        model, sampleTE,
        obsX.data(), obsY.data(), obsZ.data(), nObs,
        S, doCps,
        &rst, &fwd,
        warmStartPtr,
        aWeighting,
        ncrithyst);

    // ── pack result ───────────────────────────────────────────────────────────
    py::dict result;
    result["conv"] = converged ? 1 : 0;
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
    }
    return result;
}

PYBIND11_MODULE(gfoil_cpp, m) {
    m.def("run_forward", &run_forward_py,
          py::arg("input_dict"),
          py::arg("prev_jacobian") = py::none(),
          "Run forward aero+acoustic solver. Returns dict with CL/CD/CM/OASPL and jacobian.");
    m.def("run_AD", &run_AD_py,
          "Run AD solver given input dict and jacobian from run_forward.");
}
