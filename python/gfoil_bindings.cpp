#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <vector>

// ── Unified physics headers ────────────────────────────────────────────────
#include "codi.hpp"
using Real = codi::RealReverse;
using Tape = typename Real::Tape;

#include "real_type.hpp"
#include "data_structs.hpp"
#include "spline.hpp"
#include "main_func.hpp"
#include "get_funcs.hpp"
#include "residuals.hpp"
#include "calc_ue_m.hpp"
#include "solve_inv.hpp"
#include "build_global_sys.hpp"
#include "solve_glob.hpp"
#include "clear_RV.hpp"
#include "stagmove.hpp"
#include "coupled.hpp"
#include "init_BL.hpp"
#include "update_state.hpp"
#include "update_transition.hpp"
#include "extract_BL_TE.hpp"
#include "sound.hpp"
#include "ADfuncs.hpp"
// ─────────────────────────────────────────────────────────────────────────

namespace py = pybind11;

// ── Forward solve ──────────────────────────────────────────────────────────
py::dict run_forward(
    py::array_t<double> xcoords_np,
    py::array_t<double> ycoords_np,
    double alpha_deg,
    double Re_val, double Ma_val,
    double rho_val, double nu_val,
    double ncrit_val, double chord_val,
    double top_trans, double bot_trans,
    bool force_trans,
    double sampleTE_val,
    double obs_x, double obs_y, double obs_z, double obs_s,
    const std::string& wps_model,
    bool return_data)
{
    auto xbuf = xcoords_np.unchecked<1>();
    auto ybuf = ycoords_np.unchecked<1>();
    if (xbuf.size() != Nin || ybuf.size() != Nin)
        throw std::invalid_argument("xcoords/ycoords must each have Nin elements");

    Real inXcoords[Nin], inYcoords[Nin];
    for (int i=0;i<Nin;++i) { inXcoords[i]=xbuf(i); inYcoords[i]=ybuf(i); }

    Real alphaDeg  = alpha_deg;
    Real ReNum     = Re_val;
    Real MaNum     = Ma_val;
    Real rhoInf    = rho_val;
    Real nuInf     = nu_val;
    Real custChord = chord_val;
    Real ncrit     = ncrit_val;
    Real sampleTE  = sampleTE_val;
    const Real Ufac = 1.0, TEfac = 1.0;

    Real alpha = (alphaDeg / 180.0) * M_PI;
    Oper<Real> oper(alpha, ReNum, MaNum);
    oper.rho = rhoInf;
    Geom<Real> geom;

    Real flatCoords[2*Ncoords]={0}, inCoords[2*Nin]={0};
    for (int i=0;i<Nin;++i) {
        inCoords[colMajorIndex(0,i,2)] = inXcoords[i];
        inCoords[colMajorIndex(1,i,2)] = inYcoords[i];
    }
    make_panels<Real>(inCoords, flatCoords, Ufac, TEfac);

    Foil<Real> foil(flatCoords);
    Isolc<Real> isolc;
    Isolv<Real> isol_v;
    Param<Real> param;
    param.ncrit = ncrit;
    Wake<Real> wake;
    Vsol<Real> vsol;
    Glob<Real> glob;

    build_gamma_codi<Real>(isolc, foil, oper);
    init_thermo<Real>(oper, param, geom);
    build_wake<Real>(foil, geom, oper, isolc, wake);
    stagpoint_find<Real>(isolc, isol_v, foil, wake);
    identify_surfaces<Real>(isol_v, vsol);
    set_wake_gap<Real>(foil, isol_v, vsol);
    calc_ue_m<Real>(foil, wake, isolc, vsol);
    rebuild_ue_m<Real>(foil, wake, isol_v, vsol, false);

    Real xTransTop = Real(top_trans)  * geom.chord;
    Real xTransBot = Real(bot_trans)  * geom.chord;
    int idx_bot = vsol.Is[0].empty() ? 0 : vsol.Is[0].front();
    int idx_top = vsol.Is[1].empty() ? Ncoords-1 : vsol.Is[1].back();
    if (force_trans) {
        for (int i : vsol.Is[0]) {
            if (flatCoords[colMajorIndex(0,i,2)] >= xTransBot) { idx_bot=i; break; }
        }
        for (int i : vsol.Is[1]) {
            if (flatCoords[colMajorIndex(0,i,2)] >= xTransTop) { idx_top=i; break; }
        }
    }
    Trans<Real> tdata;
    tdata.transNode[0]=idx_bot; tdata.transNode[1]=idx_top;
    tdata.transPos[0]=xTransBot; tdata.transPos[1]=xTransTop;

    init_boundary_layer<Real>(oper, foil, param, isolc, isol_v, vsol, glob, tdata, force_trans);
    stagpoint_move<Real>(isol_v, glob, foil, wake, vsol);
    bool converged = solve_coupled<Real>(oper, foil, wake, param, vsol, isolc, isol_v, glob, tdata, force_trans);

    Post<Real> post;
    calc_force<Real>(oper, geom, param, foil, glob, post);

    Real Uinf = (ReNum * nuInf) / custChord;

    Real topsurf[7]={0}, botsurf[7]={0};
    Real xcoords_arr[Ncoords]={0};
    for (int i=0;i<Ncoords;++i) xcoords_arr[i] = flatCoords[colMajorIndex(0,i,2)];
    interpolate_at_95_both_surfaces<Real>(xcoords_arr, glob.U.data(), post.cp.data(), oper, vsol.turb.data(), param,
        topsurf, botsurf, Uinf, sampleTE, custChord);
    Real OASPL = calc_OASPL_AD<Real>(botsurf, topsurf, custChord, Uinf,
        Real(obs_x), Real(obs_y), Real(obs_z), Real(obs_s), nuInf, rhoInf, wps_model);

    if (std::isnan(OASPL.getValue()) || std::isinf(OASPL.getValue())) converged = false;

    py::dict result;
    result["converged"] = converged;
    result["CL"] = post.cl.getValue();
    result["CD"] = post.cd.getValue();
    result["CM"] = post.cm.getValue();
    result["OASPL"] = OASPL.getValue();
    result["Uinf"] = Uinf.getValue();

    if (return_data) {
        std::vector<double> inner(2*Ncoords);
        for (int i=0;i<2*Ncoords;++i) inner[i] = foil.x[i].getValue();
        std::vector<double> cp_out(Ncoords+Nwake);
        for (int i=0;i<Ncoords+Nwake;++i) cp_out[i] = post.cp[i].getValue();
        result["innerFoil"]           = inner;
        result["Cp"]                  = cp_out;
        result["thetaUpper"]          = topsurf[0].getValue();
        result["deltaStarUpper"]      = topsurf[1].getValue();
        result["tauMaxUpper"]         = topsurf[2].getValue();
        result["edgeVelocityUpper"]   = topsurf[3].getValue();
        result["dpdxUpper"]           = topsurf[4].getValue();
        result["tauWallUpper"]        = topsurf[5].getValue();
        result["delta99Upper"]        = topsurf[6].getValue();
        result["thetaLower"]          = botsurf[0].getValue();
        result["deltaStarLower"]      = botsurf[1].getValue();
        result["tauMaxLower"]         = botsurf[2].getValue();
        result["edgeVelocityLower"]   = botsurf[3].getValue();
        result["dpdxLower"]           = botsurf[4].getValue();
        result["tauWallLower"]        = botsurf[5].getValue();
        result["delta99Lower"]        = botsurf[6].getValue();
        result["stagnation"]          = std::vector<int>{isol_v.stagIndex[0], isol_v.stagIndex[1]};
    }

    return result;
}

// ── Adjoint solve ─────────────────────────────────────────────────────────
py::dict run_adjoint_solve(
    py::array_t<double> xcoords_np,
    py::array_t<double> ycoords_np,
    double alpha_deg,
    double Re_val, double Ma_val,
    double rho_val, double nu_val,
    double ncrit_val, double chord_val,
    double top_trans, double bot_trans,
    double sampleTE_val,
    double obs_x, double obs_y, double obs_z, double obs_s,
    const std::string& wps_model,
    // restart data (converged state from a prior run_forward call)
    py::array_t<double> states_np,
    py::array_t<int>    turb_np,
    py::array_t<int>    stag_np,
    py::array_t<double> rv_vals_np,
    py::array_t<int>    rv_rows_np,
    py::array_t<int>    rv_cols_np)
{
    // This mirrors srcAD/main.cpp — loads precomputed state and Jacobian,
    // then runs the adjoint computation.
    throw std::runtime_error("run_adjoint: not yet implemented in pybind11 path. "
        "Use the GFoil_AD executable with restart.json until full adjoint binding is complete.");
}

// ── Module definition ─────────────────────────────────────────────────────
PYBIND11_MODULE(gfoil_core, m) {
    m.doc() = "GFoil coupled viscous-inviscid aerofoil solver with adjoint";

    m.def("run", &run_forward,
        py::arg("xcoords"), py::arg("ycoords"),
        py::arg("alpha"), py::arg("Re"), py::arg("Ma"),
        py::arg("rho"), py::arg("nu"),
        py::arg("ncrit"), py::arg("chord"),
        py::arg("top_trans")=0.7, py::arg("bot_trans")=0.7,
        py::arg("force_trans")=false,
        py::arg("sampleTE")=0.95,
        py::arg("obs_x")=0.0, py::arg("obs_y")=3.0, py::arg("obs_z")=0.0, py::arg("obs_s")=1.0,
        py::arg("model")="Goody",
        py::arg("return_data")=false,
        "Run the forward aerofoil solve. Returns dict with CL, CD, CM, OASPL, converged.");

    m.def("run_adjoint", &run_adjoint_solve,
        py::arg("xcoords"), py::arg("ycoords"),
        py::arg("alpha"), py::arg("Re"), py::arg("Ma"),
        py::arg("rho"), py::arg("nu"),
        py::arg("ncrit"), py::arg("chord"),
        py::arg("top_trans"), py::arg("bot_trans"),
        py::arg("sampleTE"),
        py::arg("obs_x"), py::arg("obs_y"), py::arg("obs_z"), py::arg("obs_s"),
        py::arg("model"),
        py::arg("states"), py::arg("turb"), py::arg("stag"),
        py::arg("rv_vals"), py::arg("rv_rows"), py::arg("rv_cols"),
        "Run the adjoint solve given a precomputed converged state.");
}
