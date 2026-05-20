#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.h"
#include "data_structs.h"
#include "spline.hpp"
#include "main_func.h"
#include "get_funcs.hpp"
#include "extract_BL_TE.hpp"
#include "panel_funcs.hpp"
#include "calc_ue_m.hpp"
#include "solve_inv.hpp"
#include "solver_funcs.hpp"
#include "sound.hpp"
#include "restart_state.h"
#include "run_forward.h"
#include <fstream>
#include <string>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

bool runCode(
    bool fromRestart,
    const Real nCrit,
    const Real Ufac,
    const Real TEfac,
    const Real chordScaling,
    const Real (&inXcoords)[Nin],
    Real (&inYcoords)[Nin],
    Real alphad,
    Real Re,
    Real Ma,
    Real rhoInf,
    Real kinViscInf,
    const std::string model,
    const Real sampleTE,
    const Real* obsX,
    const Real* obsY,
    const Real* obsZ,
    int nObs,
    const Real S,
    const int doCps,
    RestartState* restartOut,
    ForwardResult* fwdOut,
    const RestartState* warmStart,
    int aWeighting,
    Real ncrithyst)
{
    Real alpha = (alphad / 180) * M_PI;
    Oper oper(alpha, Re, Ma);
    oper.rho = rhoInf;

    Geom geom;

    Real flattenedCoords[2 * Ncoords] = {0};
    Real inCoords[2 * Nin] = {0};
    for (int i = 0; i < Nin; ++i) {
        inCoords[colMajorIndex(0, i, 2)] = inXcoords[i];
        inCoords[colMajorIndex(1, i, 2)] = inYcoords[i];
    }
    make_panels(inCoords, flattenedCoords, Ufac, TEfac);

    Foil foil(flattenedCoords);
    Isol isol;
    Param param;
    param.ncrit     = nCrit;
    param.ncrithyst = ncrithyst;
    Wake wake;
    Vsol vsol;
    Glob glob;

    build_gamma_codi(isol, foil, oper);
    init_thermo<>(oper, param, geom);
    build_wake_impl<>(foil, geom, oper, isol, wake);
    stagpoint_find_impl<true>(isol, isol, foil, wake);
    identify_surfaces<>(isol, vsol);
    set_wake_gap<>(foil, isol, vsol);
    calc_ue_m<Real>(foil, wake, isol, vsol);
    rebuild_ue_m<>(foil, wake, isol, vsol, false);

    if (warmStart != nullptr) {
        // pybind11 path: initialise from in-memory state passed by caller
        for (int i = 0; i < RVdimension; ++i)
            glob.U[i] = warmStart->states[i];
        for (int i = 0; i < (Ncoords + Nwake); ++i)
            vsol.turb[i] = static_cast<bool>(warmStart->turb[i]);
    } else if (fromRestart) {
        // binary path: read restart.json from disk
        std::ifstream prevfile("restart.json");
        json j;
        prevfile >> j;
        for (int i = 0; i < RVdimension; ++i) {
            glob.U[i] = j["states"][i].get<double>();
        }
        for (int i = 0; i < (Ncoords + Nwake); ++i) {
            vsol.turb[i] = j["turb"][i].get<bool>();
        }
    } else {
        init_boundary_layer(oper, foil, param, isol, vsol, glob);
    }

    stagpoint_move(isol, glob, foil, wake, vsol);
    bool converged = solve_coupled(oper, foil, wake, param, vsol, isol, glob, restartOut);
    Post post;
    calc_force<>(oper, geom, param, foil, glob, post);

    Real Uinf = (Re * kinViscInf) / (chordScaling);
    Real tauWall[Ncoords];
    if (doCps) {
        Real cf_U[4] = {0};
        for (int i = 0; i < Ncoords; ++i) {
            tauWall[i] = get_cf(
                glob.U[colMajorIndex(0, i, 4)],
                glob.U[colMajorIndex(1, i, 4)],
                glob.U[colMajorIndex(2, i, 4)],
                glob.U[colMajorIndex(3, i, 4)],
                vsol.turb[i], false, param, cf_U);
            tauWall[i] *= (oper.rho * (glob.U[colMajorIndex(3, i, 4)] * Uinf *
                                       glob.U[colMajorIndex(3, i, 4)] * Uinf)) / 2;
        }
    }

    Real topsurf[7], botsurf[7];
    Real xcoords[Ncoords] = {0};
    Real ycoords[Ncoords] = {0};
    for (int i = 0; i < Ncoords; ++i) {
        xcoords[i] = flattenedCoords[colMajorIndex(0, i, 2)];
        ycoords[i] = flattenedCoords[colMajorIndex(1, i, 2)];
    }

    interpolate_at_95_both_surfaces(xcoords, glob.U, post.cp, oper, vsol.turb,
                                    param, topsurf, botsurf, Uinf, sampleTE, chordScaling);
    Real OASPL = calc_OASPL<Real, true>(botsurf, topsurf, chordScaling, Uinf,
                                         obsX, obsY, obsZ, nObs, S, kinViscInf, rhoInf, model, doCps, aWeighting);

    if (std::isnan(OASPL) || std::isinf(OASPL))
        converged = false;

    if (fwdOut != nullptr) {
        // Pybind11 path: fill result struct, skip file writes
        fwdOut->converged = converged;
        if (converged) {
            fwdOut->CL    = post.cl.getValue();
            fwdOut->CD    = post.cd.getValue();
            fwdOut->CM    = post.cm.getValue();
            fwdOut->OASPL = OASPL.getValue();
        }
        return converged;
    }

    // Standalone binary path: write out.json
    if (converged) {
        json out;
        out["conv"] = 1;
        out["aerofoilChord"]      = chordScaling.getValue();
        out["freestreamVelocity"] = Uinf.getValue();
        out["samplingLoc"]        = sampleTE.getValue();
        out["CL"]    = post.cl.getValue();
        out["CD"]    = post.cd.getValue();
        out["CM"]    = post.cm.getValue();
        out["OASPL"] = OASPL.getValue();

        if (doCps) {
            Real botTransX = geom.chord;
            for (int i = 0; i < isol.stagIndex[0]; ++i) {
                if (vsol.turb[isol.stagIndex[0] - i]) {
                    botTransX = foil.x[colMajorIndex(0, isol.stagIndex[0] - i, 2)];
                    break;
                }
            }
            Real topTransX = geom.chord;
            for (int i = 0; i < 200 - isol.stagIndex[1]; ++i) {
                if (vsol.turb[isol.stagIndex[1] + i]) {
                    topTransX = foil.x[colMajorIndex(0, isol.stagIndex[1] + i, 2)];
                    break;
                }
            }

            double inner[2 * Ncoords];
            double cps[Ncoords];
            double tauWallOut[Ncoords];
            for (int i = 0; i < 2 * Ncoords; ++i) inner[i] = foil.x[i].getValue();
            for (int i = 0; i < Ncoords; ++i) {
                cps[i]        = post.cp[i].getValue();
                tauWallOut[i] = tauWall[i].getValue();
            }
            out["innerFoil"]   = inner;
            out["Cp"]          = cps;
            out["tauWall"]     = tauWallOut;
            out["stagnation"]  = isol.stagIndex;
            out["topTransX"]   = topTransX.getValue();
            out["botTransX"]   = botTransX.getValue();

            out["thetaUpper"]        = topsurf[0].getValue();
            out["deltaStarUpper"]    = topsurf[1].getValue();
            out["tauMaxUpper"]       = topsurf[2].getValue();
            out["edgeVelocityUpper"] = topsurf[3].getValue();
            out["dpdxUpper"]         = topsurf[4].getValue();
            out["tauWallUpper"]      = topsurf[5].getValue();
            out["delta99Upper"]      = topsurf[6].getValue();

            out["thetaLower"]        = botsurf[0].getValue();
            out["deltaStarLower"]    = botsurf[1].getValue();
            out["tauMaxLower"]       = botsurf[2].getValue();
            out["edgeVelocityLower"] = botsurf[3].getValue();
            out["dpdxLower"]         = botsurf[4].getValue();
            out["tauWallLower"]      = botsurf[5].getValue();
            out["delta99Lower"]      = botsurf[6].getValue();
        }

        std::ofstream outFile("out.json");
        outFile << out.dump(4);
    } else {
        json out;
        out["conv"] = 0;
        std::ofstream outFile("out.json");
        outFile << out.dump(4);
    }

    return converged;
}
