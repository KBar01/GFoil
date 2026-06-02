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
#include <string>
#include <memory>

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
    RestartState* restartOut,
    ForwardResult* fwdOut,
    const RestartState* warmStart,
    int aWeighting,
    Real ncrithyst,
    bool verbose,
    double f_min,
    double f_max,
    double xft_lower,
    double xft_upper,
    double rtol)
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
    auto isolPtr = std::make_unique<Isol>();
    auto vsolPtr = std::make_unique<Vsol>();
    auto globPtr = std::make_unique<Glob>();
    Isol& isol = *isolPtr;
    Vsol& vsol = *vsolPtr;
    Glob& glob = *globPtr;
    Param param;
    param.rtol      = rtol;       // RMS convergence tolerance (forward-only knob)
    param.ncrit     = nCrit;
    param.ncrithyst = ncrithyst;
    param.xft_xc[0] = xft_lower;  // lower surface (vsol.Is[0])
    param.xft_xc[1] = xft_upper;  // upper surface (vsol.Is[1])
    Wake wake;

    build_gamma_codi(isol, foil, oper);
    init_thermo<>(oper, param, geom);
    build_wake_impl<>(foil, geom, oper, isol, wake);
    stagpoint_find_impl<true>(isol, isol, foil, wake);
    identify_surfaces<>(isol, vsol);
    set_wake_gap<>(foil, isol, vsol);
    calc_ue_m<Real>(foil, wake, isol, vsol);
    rebuild_ue_m<>(foil, wake, isol, vsol, false);

    if (warmStart != nullptr) {
        for (int i = 0; i < RVdimension; ++i)
            glob.U[i] = warmStart->states[i];
        for (int i = 0; i < (Ncoords + Nwake); ++i)
            vsol.turb[i] = static_cast<bool>(warmStart->turb[i]);
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

    Real topsurf[7], botsurf[7];
    Real xcoords[Ncoords] = {0};
    Real ycoords[Ncoords] = {0};
    for (int i = 0; i < Ncoords; ++i) {
        xcoords[i] = flattenedCoords[colMajorIndex(0, i, 2)];
        ycoords[i] = flattenedCoords[colMajorIndex(1, i, 2)];
    }

    interpolate_at_95_both_surfaces(xcoords, glob.U, post.cp, oper, vsol.turb,
                                    param, topsurf, botsurf, Uinf, sampleTE, chordScaling);
    Real OASPL = calc_OASPL<Real>(botsurf, topsurf, chordScaling, Uinf,
                                   obsX, obsY, obsZ, nObs, S, kinViscInf, rhoInf, model,
                                   f_min, f_max, aWeighting, alpha);

    // The aerodynamic solve can converge while the downstream acoustic model
    // (Amiet/WPS) yields a non-finite OASPL — common at low Re where the BL
    // edge quantities feeding the noise model degenerate.  Distinguish this
    // from a genuine convergence failure instead of reporting a silent blank
    // failure_mode that masquerades as non-convergence.  (Noise code untouched.)
    const bool aero_converged = converged;
    const bool acoustic_nan   = (std::isnan(OASPL) || std::isinf(OASPL));
    if (acoustic_nan) converged = false;

    if (fwdOut != nullptr) {
        fwdOut->converged          = converged;
        fwdOut->failure_mode       = converged ? ""
            : (aero_converged && acoustic_nan ? "acoustic_nan" : failure_mode);
        fwdOut->newton_iterations  = glob.convergenceIteration;
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
                    Real log_fmin = std::log10(Real(f_min));
                    Real log_fmax = std::log10(Real(f_max));
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

                    const double pref2   = (20e-6) * (20e-6);
                    const double two_pi  = 2.0 * M_PI;
                    for (int i = 0; i < NS; ++i) {
                        double wps_u_hz = WPS_U_arr[i].getValue() * two_pi;
                        double wps_l_hz = WPS_L_arr[i].getValue() * two_pi;
                        fwdOut->WPS_upper[i] = (wps_u_hz > 0.0) ? 10.0 * std::log10(wps_u_hz / pref2) : -200.0;
                        fwdOut->WPS_lower[i] = (wps_l_hz > 0.0) ? 10.0 * std::log10(wps_l_hz / pref2) : -200.0;
                    }

                    // Transform to TE-local Amiet frame (consistent with calc_OASPL)
                    const Real cos_a     = std::cos(alpha);
                    const Real sin_a     = std::sin(alpha);
                    const Real te_offset = static_cast<Real>(0.75) * chordScaling;
                    for (int iObs = 0; iObs < nObs; ++iObs) {
                        Real ff[Nsound];
                        for (int i = 0; i < NS; ++i) ff[i] = 0.0;
                        Real x_loc = obsX[iObs] * cos_a - obsZ[iObs] * sin_a - te_offset;
                        Real y_loc = obsY[iObs];
                        Real z_loc = obsX[iObs] * sin_a + obsZ[iObs] * cos_a;
                        TE_noise_outer<Real>(Uinf / 340.0, Uinf,
                                             x_loc, y_loc, z_loc,
                                             chordScaling / 2.0, chordScaling,
                                             S, Real(340.0), omArr,
                                             edgeVel_bot, edgeVel_top,
                                             WPS_L_arr, WPS_U_arr, ff);
                        for (int i = 0; i < NS; ++i) {
                            double ff_hz = ff[i].getValue() * two_pi;
                            fwdOut->FF_spectra[iObs * NS + i] = (ff_hz > 0.0) ? 10.0 * std::log10(ff_hz / pref2) : -200.0;
                        }
                    }
                }
            } // end verbose block
        }
        return converged;
    }

    return converged;
}
