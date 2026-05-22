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
    Real ncrithyst,
    bool verbose)
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
    std::string failure_mode;
    bool converged = solve_coupled(oper, foil, wake, param, vsol, isol, glob, restartOut,
                                   (fwdOut != nullptr) ? &failure_mode : nullptr);
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
        fwdOut->converged     = converged;
        fwdOut->failure_mode  = converged ? "" : failure_mode;
        if (converged) {
            fwdOut->CL    = post.cl.getValue();
            fwdOut->CD    = post.cd.getValue();
            fwdOut->CM    = post.cm.getValue();
            fwdOut->OASPL = OASPL.getValue();

            if (verbose) {
                const int N = Ncoords;
                fwdOut->innerFoilX.resize(N);
                fwdOut->innerFoilY.resize(N);
                fwdOut->Cp.resize(N);
                fwdOut->delta_star.resize(N);
                fwdOut->theta.resize(N);
                fwdOut->tau_wall.resize(N);
                fwdOut->tau_max.resize(N);
                fwdOut->Ue.resize(N);
                fwdOut->dpdx.resize(N);
                fwdOut->is_turb.resize(N);

                Real cf_U[4] = {0};
                double chord = chordScaling.getValue();

                for (int i = 0; i < N; ++i) {
                    fwdOut->innerFoilX[i] = xcoords[i].getValue();
                    fwdOut->innerFoilY[i] = ycoords[i].getValue();
                    fwdOut->Cp[i]         = post.cp[i].getValue();
                    fwdOut->is_turb[i]    = vsol.turb[i];

                    double th        = glob.U[colMajorIndex(0, i, 4)].getValue();
                    double ds        = glob.U[colMajorIndex(1, i, 4)].getValue();
                    double ctau_sqrt = glob.U[colMajorIndex(2, i, 4)].getValue();

                    fwdOut->theta[i]      = th * chord;
                    fwdOut->delta_star[i] = ds * chord;

                    Real uk_ignore;
                    double Ue_phys = get_uk(glob.U[colMajorIndex(3, i, 4)], param, uk_ignore).getValue()
                                     * Uinf.getValue();
                    fwdOut->Ue[i] = Ue_phys;

                    if (vsol.turb[i]) {
                        double Ctau = ctau_sqrt * ctau_sqrt;
                        fwdOut->tau_max[i] = Ctau * oper.rho.getValue() * Ue_phys * Ue_phys;
                    } else {
                        fwdOut->tau_max[i] = 0.0;
                    }

                    double Cf = get_cf(
                        glob.U[colMajorIndex(0, i, 4)],
                        glob.U[colMajorIndex(1, i, 4)],
                        glob.U[colMajorIndex(2, i, 4)],
                        glob.U[colMajorIndex(3, i, 4)],
                        vsol.turb[i], false, param, cf_U).getValue();
                    fwdOut->tau_wall[i] = (Cf / 2.0) * oper.rho.getValue() * Ue_phys * Ue_phys;

                    int im = (i > 0)   ? i - 1 : i;
                    int ip = (i < N-1) ? i + 1 : i;
                    double dx = xcoords[ip].getValue() - xcoords[im].getValue();
                    double dCp_dx = (dx != 0.0)
                        ? (post.cp[ip].getValue() - post.cp[im].getValue()) / dx
                        : 0.0;
                    fwdOut->dpdx[i] = dCp_dx / chord * 0.5 * oper.rho.getValue()
                                      * Uinf.getValue() * Uinf.getValue();
                }

                // Transition x-locations (same coordinate frame as foil.x)
                fwdOut->botTransX = chordScaling.getValue();
                for (int i = 0; i < isol.stagIndex[0]; ++i) {
                    if (vsol.turb[isol.stagIndex[0] - i]) {
                        fwdOut->botTransX = foil.x[colMajorIndex(0, isol.stagIndex[0] - i, 2)].getValue();
                        break;
                    }
                }
                fwdOut->topTransX = chordScaling.getValue();
                for (int i = 0; i < 200 - isol.stagIndex[1]; ++i) {
                    if (vsol.turb[isol.stagIndex[1] + i]) {
                        fwdOut->topTransX = foil.x[colMajorIndex(0, isol.stagIndex[1] + i, 2)].getValue();
                        break;
                    }
                }

                // TE BL sampling states
                fwdOut->BL_top.resize(7);
                fwdOut->BL_bot.resize(7);
                for (int i = 0; i < 7; ++i) {
                    fwdOut->BL_top[i] = topsurf[i].getValue();
                    fwdOut->BL_bot[i] = botsurf[i].getValue();
                }

                // Acoustic spectra — rebuild freq/WPS/FF from the same BL states
                {
                    const int NS = Nsound;
                    fwdOut->freq_Hz.resize(NS);
                    fwdOut->WPS_upper.resize(NS, 0.0);
                    fwdOut->WPS_lower.resize(NS, 0.0);
                    fwdOut->nObs = nObs;
                    fwdOut->FF_spectra.resize(nObs * NS, 0.0);

                    Real omArr[Nsound];
                    Real Freq_arr[Nsound];
                    Real log_fmin = std::log10(Real(200.0));
                    Real log_fmax = std::log10(Real(20000.0));
                    for (int i = 0; i < NS; ++i) {
                        Real frac  = Real(i) / Real(NS - 1);
                        Real logf  = log_fmin + frac * (log_fmax - log_fmin);
                        Freq_arr[i] = std::pow(Real(10.0), logf);
                        omArr[i]    = 2.0 * M_PI * Freq_arr[i];
                        fwdOut->freq_Hz[i] = Freq_arr[i].getValue();
                    }

                    Real WPS_U_arr[Nsound];
                    Real WPS_L_arr[Nsound];
                    for (int i = 0; i < NS; ++i) { WPS_U_arr[i] = 0.0; WPS_L_arr[i] = 0.0; }

                    // Upper surface WPS (mirrors calc_OASPL logic)
                    Real tauWall_top = topsurf[5];
                    if (tauWall_top < 0.0) tauWall_top *= -1.0;
                    Real edgeVel_top = topsurf[3];
                    if (topsurf[2] > 0.0) {
                        calc_WPS<Real>(model, topsurf[0], topsurf[1], topsurf[6],
                                       tauWall_top, topsurf[2], edgeVel_top, topsurf[4],
                                       omArr, kinViscInf, Uinf,
                                       obsX[0], obsY[0], obsZ[0], S, rhoInf, 1, WPS_U_arr);
                    } else {
                        edgeVel_top = Uinf;
                    }

                    // Lower surface WPS
                    Real tauWall_bot = botsurf[5];
                    if (tauWall_bot < 0.0) tauWall_bot *= -1.0;
                    Real edgeVel_bot = botsurf[3];
                    if (botsurf[2] > 0.0) {
                        calc_WPS<Real>(model, botsurf[0], botsurf[1], botsurf[6],
                                       tauWall_bot, botsurf[2], edgeVel_bot, botsurf[4],
                                       omArr, kinViscInf, Uinf,
                                       obsX[0], obsY[0], obsZ[0], S, rhoInf, 0, WPS_L_arr);
                    } else {
                        edgeVel_bot = Uinf;
                    }

                    for (int i = 0; i < NS; ++i) {
                        fwdOut->WPS_upper[i] = WPS_U_arr[i].getValue();
                        fwdOut->WPS_lower[i] = WPS_L_arr[i].getValue();
                    }

                    // Far-field spectrum per observer
                    for (int iObs = 0; iObs < nObs; ++iObs) {
                        Real ff[Nsound];
                        for (int i = 0; i < NS; ++i) ff[i] = 0.0;
                        TE_noise_outer<Real>(Uinf / 340.0, Uinf,
                                             obsX[iObs], obsY[iObs], obsZ[iObs],
                                             chordScaling / 2.0, chordScaling,
                                             S, Real(340.0), omArr,
                                             edgeVel_bot, edgeVel_top,
                                             WPS_L_arr, WPS_U_arr, ff);
                        for (int i = 0; i < NS; ++i)
                            fwdOut->FF_spectra[iObs * NS + i] = ff[i].getValue();
                    }
                }
            } // end verbose block
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
