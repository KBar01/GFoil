#pragma once

// Shared acoustic post-processing header.
//
// calc_WPS  : WPS model dispatch (template<typename Real>)
//             Includes all five models: roz, goo, lee, kam, tno.
//
// calc_OASPL: Compute overall SPL from boundary-layer states.
//             template<typename Real>

#include <cmath>
#include "WPSmodels.hpp"
#include "newAmiet.hpp"

// ── WPS model dispatch ────────────────────────────────────────────────────────

template<typename Real>
void calc_WPS(const std::string& model,
              const Real theta_in, const Real deltaStar_in, const Real delta_in,
              const Real tauW_in, const Real tauMax_in,
              const Real edgeVel, const Real dpdx,
              const Real (&omega)[Nsound],
              Real nu, const Real Uinf,
              const Real X, const Real Y, const Real Z,
              const Real S, const Real rho,
              const int isSuction,
              Real (&WPS)[Nsound])
{
    // ── Physically-motivated input floors (low-Re / near-separation guard) ──
    //
    // The empirical WPS models degenerate when the trailing-edge BL is mostly
    // laminar with a vanishingly thin, near-separated turbulent layer (low Re):
    // tauWall -> 0 sends u_tau -> 0 and the Clauser parameter
    // beta_c = (theta/tauWall)*dpdx -> O(100-1000), which drives the Rozenberg
    // amplitude exponent A1 = 3.7 + 1.5*beta_c past the double overflow of
    // pow(base, A1) (confirmed by the GFOIL_DEBUG trace).  deltaStar -> theta
    // likewise sends Delta=delta/deltaStar and SS=Ue/(tauMax^2*deltaStar) out of
    // range.  Floor the *physical inputs* at the thinnest resolvable attached
    // turbulent layer rather than clamping the output spectrum (which would
    // silently fabricate a noise level).  Floors are smooth std::max and, by
    // construction, never bind for a normal attached TE BL (so moderate/high-Re
    // results are bit-identical) -- verified under GFOIL_DEBUG.
    constexpr double Cf_min   = 1e-4;   // min skin-friction coefficient (well below
                                        //   any attached turbulent Cf ~ 1e-3..5e-3)
    constexpr double beta_max = 50.0;   // Clauser-parameter ceiling: edge of the
                                        //   empirical APG calibration; beyond it the
                                        //   models extrapolate meaninglessly
    constexpr double H_min    = 1.05;   // a turbulent BL has deltaStar > theta

    const Real Ue        = std::max(edgeVel,    Real(1e-6));
    const Real theta     = std::max(theta_in,   Real(1e-12));
    const Real deltaStar = std::max(deltaStar_in, Real(H_min) * theta);
    const Real delta     = std::max(delta_in,   deltaStar);
    const Real q         = Real(0.5) * rho * Ue * Ue;          // dynamic pressure
    // Min wall shear: the larger of a minimum-Cf floor and the value that keeps
    // the Clauser parameter within the model's validity ceiling.
    Real tauW            = std::max(tauW_in, Real(Cf_min) * q);
    tauW                 = std::max(tauW, theta * std::abs(dpdx) / Real(beta_max));
    const Real tauMax    = std::max(tauMax_in, tauW);          // tauMax >= tauWall

    Real useTauW = tauW;
    if (tauW > tauMax) { useTauW = tauMax; }

    if      (model == "roz") { calc_WPS_Rozenburg<Real>(theta,deltaStar,delta,useTauW,tauMax,Ue,dpdx,omega,rho,nu,WPS); }
    else if (model == "goo") { calc_WPS_Goody<Real>    (theta,deltaStar,delta,useTauW,tauMax,Ue,dpdx,omega,rho,nu,Uinf,WPS); }
    else if (model == "lee") { calc_WPS_Lee<Real>       (theta,deltaStar,delta,useTauW,tauMax,Ue,dpdx,omega,rho,nu,WPS); }
    else if (model == "kam") { calc_WPS_Kamruzzaman<Real>(theta,deltaStar,useTauW,Ue,dpdx,omega,rho,nu,WPS); }
    else if (model == "tno") { calc_WPS_TNO<Real>       (delta,useTauW,Ue,omega,rho,nu,isSuction,WPS); }
}

// ── OASPL ─────────────────────────────────────────────────────────────────────
// Overall sound pressure level integrated from PSD over Nsound frequencies
// (200-20000 Hz, log-spaced). OASPL = 10*log10(integral(PSD*domega) / pref^2)
// where pref = 20e-6 Pa.

template<typename Real>
Real calc_OASPL(const Real* botStates, const Real* topStates,
                const Real chordScale, const Real Uinf,
                const Real* obsX, const Real* obsY, const Real* obsZ,
                int nObs,
                const Real S, const Real nu, const Real rho,
                const std::string& model,
                double f_min = 200.0, double f_max = 20000.0,
                const int aWeighting = 0,
                Real alpha_rad = Real(0.0))
{
    Real omega[Nsound];
    Real Freq[Nsound];

    Real log_fmin = std::log10(static_cast<Real>(f_min));
    Real log_fmax = std::log10(static_cast<Real>(f_max));
    for (int i = 0; i < Nsound; ++i) {
        Real frac = static_cast<Real>(i) / static_cast<Real>(Nsound - 1);
        Real logf = log_fmin + frac * (log_fmax - log_fmin);
        Freq[i]  = std::pow(static_cast<Real>(10.0), logf);
        omega[i] = 2.0 * M_PI * Freq[i];
    }

    Real WPSUpper[Nsound] = {0};
    Real WPSLower[Nsound] = {0};

    // ── top surface (observer-independent) ───────────────────────────────────
    Real theta    = topStates[0];
    Real deltaS   = topStates[1];
    Real tauMax   = topStates[2];
    Real edgeVel_top = topStates[3];
    Real dpdx     = topStates[4];
    Real tauWall  = topStates[5];
    Real delta    = topStates[6];

    if (tauWall < 0.0) { tauWall *= -1.0; }

    if (tauMax > 0.0) {
        calc_WPS<Real>(model, theta, deltaS, delta, tauWall, tauMax,
                       edgeVel_top, dpdx, omega, nu, Uinf, obsX[0], obsY[0], obsZ[0], S, rho, 1, WPSUpper);
    } else {
        edgeVel_top = Uinf;
    }

    // ── bottom surface (observer-independent) ─────────────────────────────────
    theta    = botStates[0];
    deltaS   = botStates[1];
    tauMax   = botStates[2];
    Real edgeVel_bot = botStates[3];
    dpdx     = botStates[4];
    tauWall  = botStates[5];
    delta    = botStates[6];

    if (tauWall < 0.0) { tauWall *= -1.0; }

    if (tauMax > 0.0) {
        calc_WPS<Real>(model, theta, deltaS, delta, tauWall, tauMax,
                       edgeVel_bot, dpdx, omega, nu, Uinf, obsX[0], obsY[0], obsZ[0], S, rho, 0, WPSLower);
    } else {
        edgeVel_bot = Uinf;
    }

    // ── per-observer loop: far-field PSD → integrate → power average ─────────
    Real pref2 = (20e-6) * (20e-6);
    Real powerSum = 0.0;

    // Global → TE-local chord-aligned coordinate transformation.
    // Observer coords are given in global frame (origin at 1/4-chord,
    // x=freestream, z=up). TE-local frame has origin at TE, x1 along
    // chord (positive downstream), z1 normal to chord (suction-side positive).
    Real cos_a     = std::cos(alpha_rad);
    Real sin_a     = std::sin(alpha_rad);
    Real te_offset = static_cast<Real>(0.75) * chordScale;

    for (int iObs = 0; iObs < nObs; ++iObs) {
        // Transform observer from global to TE-local Amiet frame.
        Real x_loc = obsX[iObs] * cos_a - obsZ[iObs] * sin_a - te_offset;
        Real y_loc = obsY[iObs];
        Real z_loc = obsX[iObs] * sin_a + obsZ[iObs] * cos_a;

        Real farfieldSpectra[Nsound];
        Real c = Uinf / 340.0;
        TE_noise_outer<Real>(c, Uinf, x_loc, y_loc, z_loc,
                             chordScale / 2.0, chordScale,
                             S, 340.0, omega,
                             edgeVel_bot, edgeVel_top,
                             WPSLower, WPSUpper, farfieldSpectra);

        if (aWeighting) {
            for (int i = 0; i < Nsound; ++i) {
                Real f2 = Freq[i] * Freq[i];
                Real RA = (static_cast<Real>(12194.0 * 12194.0) * f2 * f2)
                        / ( (f2 + static_cast<Real>(20.6  * 20.6))
                          * std::sqrt((f2 + static_cast<Real>(107.7 * 107.7))
                                     * (f2 + static_cast<Real>(737.9 * 737.9)))
                          * (f2 + static_cast<Real>(12194.0 * 12194.0)) );
                farfieldSpectra[i] *= RA * RA;
            }
        }

        Real integral = 0.0;
        for (int i = 0; i < Nsound - 1; ++i) {
            Real df = Freq[i+1] - Freq[i];
            integral += 0.5 * (farfieldSpectra[i] + farfieldSpectra[i+1]) * 2.0 * M_PI * df;
        }

        // Acoustic floor: when both surfaces are fully laminar at the TE there
        // is no TBL-TE noise source, so the integrated mean-square pressure is
        // exactly zero and 10*log10(0) = -inf.  That is the physically correct
        // "no modelled source" limit; report a finite, clearly-silent level
        // (~-300 dB) instead of -inf (which run_forward would treat as
        // acoustic_nan).  Branch on the passive value so any source-carrying
        // case takes the *exact* original arithmetic (bit-identical tape and
        // adjoint — std::max here perturbs the OASPL gradient at ~1e-7); the
        // floor branch is a constant, correctly contributing zero sensitivity
        // (no source ⇒ no acoustic shape-sensitivity).
        Real OASPL_i;
        if ((integral / pref2).getValue() > 1e-30)
            OASPL_i = 10.0 * std::log10(integral / pref2);
        else
            OASPL_i = Real(10.0 * std::log10(1e-30));
        powerSum += std::pow(static_cast<Real>(10.0), OASPL_i / 10.0);

        if (std::getenv("GFOIL_DEBUG")) {
            auto nf = [](double v){ return !std::isfinite(v); };
            bool wU=false,wL=false,ff=false;
            for (int i=0;i<Nsound;++i){ if(nf(WPSUpper[i].getValue()))wU=true;
                if(nf(WPSLower[i].getValue()))wL=true; if(nf(farfieldSpectra[i].getValue()))ff=true; }
            std::cerr << "[OASPL obs="<<iObs<<"] Ue_bot="<<edgeVel_bot.getValue()
                      <<" Ue_top="<<edgeVel_top.getValue()
                      <<" WPSU_nan="<<wU<<" WPSL_nan="<<wL<<" FF_nan="<<ff
                      <<" integral="<<integral.getValue()
                      <<" OASPL_i="<<OASPL_i.getValue()<<"\n";
        }
    }

    Real mean_power = powerSum / static_cast<Real>(nObs);
    if (mean_power.getValue() > 1e-30)
        return 10.0 * std::log10(mean_power);
    return Real(10.0 * std::log10(1e-30));   // all observers silent: finite floor
}
