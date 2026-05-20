#pragma once

// Shared acoustic post-processing header.
//
// calc_WPS  : WPS model dispatch (template<typename Real>)
//             Includes all five models: roz, goo, lee, kam, tno.
//
// calc_OASPL: Compute overall SPL from boundary-layer states.
//             template<typename Real, bool WriteJSON = false>
//             When WriteJSON=true  the JSON output block is compiled in
//               (used by the forward solver, instantiated in src/sound.cpp).
//             When WriteJSON=false the block is compiled out entirely
//               (used by the AD solver — zero overhead).
//             The WPSjson runtime flag is only meaningful when WriteJSON=true;
//             the if constexpr ensures it is never evaluated when WriteJSON=false.

#include <cmath>
#include <fstream>
#include <sstream>
#include "WPSmodels.hpp"
#include "newAmiet.hpp"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

// ── WPS model dispatch ────────────────────────────────────────────────────────

template<typename Real>
void calc_WPS(const std::string& model,
              const Real theta, const Real deltaStar, const Real delta,
              const Real tauW, const Real tauMax,
              const Real edgeVel, const Real dpdx,
              const Real (&omega)[Nsound],
              Real nu, const Real Uinf,
              const Real X, const Real Y, const Real Z,
              const Real S, const Real rho,
              const int isSuction,
              Real (&WPS)[Nsound])
{
    Real useTauW = tauW;
    if (tauW > tauMax) { useTauW = tauMax; }

    if      (model == "roz") { calc_WPS_Rozenburg<Real>(theta,deltaStar,delta,useTauW,tauMax,edgeVel,dpdx,omega,rho,nu,WPS); }
    else if (model == "goo") { calc_WPS_Goody<Real>    (theta,deltaStar,delta,useTauW,tauMax,edgeVel,dpdx,omega,rho,nu,Uinf,WPS); }
    else if (model == "lee") { calc_WPS_Lee<Real>       (theta,deltaStar,delta,useTauW,tauMax,edgeVel,dpdx,omega,rho,nu,WPS); }
    else if (model == "kam") { calc_WPS_Kamruzzaman<Real>(theta,deltaStar,useTauW,edgeVel,dpdx,omega,rho,nu,WPS); }
    else if (model == "tno") { calc_WPS_TNO<Real>       (delta,useTauW,edgeVel,omega,rho,nu,isSuction,WPS); }
}

// ── OASPL ─────────────────────────────────────────────────────────────────────
// Overall sound pressure level integrated from PSD over Nsound frequencies
// (200-20000 Hz, log-spaced). OASPL = 10*log10(integral(PSD*domega) / pref^2)
// where pref = 20e-6 Pa.

template<typename Real, bool WriteJSON = false>
Real calc_OASPL(const Real* botStates, const Real* topStates,
                const Real chordScale, const Real Uinf,
                const Real* obsX, const Real* obsY, const Real* obsZ,
                int nObs,
                const Real S, const Real nu, const Real rho,
                const std::string& model,
                const int WPSjson = 0,
                const int aWeighting = 0)
{
    const Real f_min = 200.0;
    const Real f_max = 20000.0;

    Real omega[Nsound];
    Real Freq[Nsound];

    Real log_fmin = std::log10(f_min);
    Real log_fmax = std::log10(f_max);
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
    Real farfieldSpectra0[Nsound];  // first observer's spectra, cached for JSON

    for (int iObs = 0; iObs < nObs; ++iObs) {
        Real farfieldSpectra[Nsound];
        Real c = Uinf / 340.0;
        TE_noise_outer<Real>(c, Uinf, obsX[iObs], obsY[iObs], obsZ[iObs],
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

        if (iObs == 0) {
            for (int i = 0; i < Nsound; ++i) { farfieldSpectra0[i] = farfieldSpectra[i]; }
        }

        Real integral = 0.0;
        for (int i = 0; i < Nsound - 1; ++i) {
            Real df = Freq[i+1] - Freq[i];
            integral += 0.5 * (farfieldSpectra[i] + farfieldSpectra[i+1]) * 2.0 * M_PI * df;
        }

        Real OASPL_i = 10.0 * std::log10(integral / pref2);
        powerSum += std::pow(static_cast<Real>(10.0), OASPL_i / 10.0);
    }

    Real OASPL = 10.0 * std::log10(powerSum / static_cast<Real>(nObs));

    // ── optional JSON output (compiled out when WriteJSON=false) ──────────────
    if constexpr (WriteJSON) {
        if (WPSjson == 1) {
            json j;
            double prefSqrd = pref2.getValue();
            std::vector<double> freq_d(Nsound);
            std::vector<double> wpsupper_d(Nsound);
            std::vector<double> wpslower_d(Nsound);
            std::vector<double> spectra_d(Nsound);
            for (int i = 0; i < Nsound; ++i) {
                freq_d[i]     = Freq[i].getValue();
                wpsupper_d[i] = WPSUpper[i].getValue() / prefSqrd;
                wpslower_d[i] = WPSLower[i].getValue() / prefSqrd;
                spectra_d[i]  = farfieldSpectra0[i].getValue() / prefSqrd;
            }
            j["frequency_Hz"]        = freq_d;
            j["WPS_upper/prefSqrd"]  = wpsupper_d;
            j["WPS_lower/prefSqrd"]  = wpslower_d;
            j["FF_spectra/prefSqrd"] = spectra_d;
            j["OASPL_dB"]            = OASPL.getValue();

            std::ofstream file("WPS.json");
            file << j.dump(4);
            file.close();
        }
    }

    return OASPL;
}
