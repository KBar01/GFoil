#pragma once

// noise_run — acoustics-only forward investigation entry point.
//
// Runs the WPS models (+ Amiet far-field model) directly from supplied
// trailing-edge boundary-layer states (or a custom wall-pressure spectrum),
// with NO aerodynamic solve. Uses the same Real = codi::RealReverse type as the
// forward path so the existing calc_WPS_*/Amiet code (which needs a CoDiPack
// active type for errFunc) is reused without a double retype — but this function
// is NEVER AD'd. No gradient is ever taken through it.
//
// Length-generic: accepts a frequency array of ANY length / ANY spacing via the
// *_vec acoustic overloads (calc_WPS_vec / TE_noise_outer_vec). The fixed-Nsound
// templates on the AD-critical forward path are untouched.
//
// All returned spectra are RAW LINEAR Pa^2/omega (no dB, no integration, no OASPL).

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include "sound.hpp"

// Plain (non-template) result holder — values are extracted to double.
struct NoiseRunResult {
    std::vector<double> freqs_Hz;        // N (echo of input)
    std::vector<double> WPS_upper;       // N      Pa^2/omega (0 if surface skipped)
    std::vector<double> WPS_lower;       // N      Pa^2/omega (0 if surface skipped)
    std::vector<double> FF_spectra;      // nObs*N Pa^2/omega, row-major (obs, freq)
    std::vector<double> obsXYZ_TElocal;  // nObs*3 observer coords in TE-local frame
};

// BL states are passed as [upper, lower] pairs (index 0 = upper, 1 = lower),
// matching the calc_WPS contract ordering used by topsurf/botsurf:
//   theta, delta_star, tau_max, Ue, dpdx, tau_wall, delta99.
template<typename Real>
NoiseRunResult noise_run_cpp(
    double alphaDeg,
    double Re, double rho, double nu, double Ma, double chord,
    const std::vector<double>& obsX,
    const std::vector<double>& obsY,
    const std::vector<double>& obsZ,
    double span,
    const std::array<double,2>& theta,
    const std::array<double,2>& deltaStar,
    const std::array<double,2>& tauMax,
    const std::array<double,2>& Ue,
    const std::array<double,2>& dpdx,
    const std::array<double,2>& tauWall,
    const std::array<double,2>& delta99,
    const std::vector<double>& freqs_Hz,
    const std::string& model,
    bool has_custom,
    const std::vector<double>& custom_upper,   // length N (or empty)
    const std::vector<double>& custom_lower)   // length N (or empty)
{
    // ── Tape hygiene ──────────────────────────────────────────────────────────
    // errFunc (in the Amiet path) pushes CoDiPack statements via
    // StatementPushHelper. We never evaluate gradients here, but rewind the tape
    // at the start of every call to prevent unbounded tape memory growth across
    // repeated in-process noise_run calls. Mirrors the reset pattern in
    // gfoil_ad_bindings.cpp's run_AD_py.
    Real::getTape().reset();

    const std::size_t N    = freqs_Hz.size();
    const std::size_t nObs = obsX.size();

    NoiseRunResult out;
    out.freqs_Hz       = freqs_Hz;
    out.WPS_upper.assign(N, 0.0);
    out.WPS_lower.assign(N, 0.0);
    out.FF_spectra.assign(nObs * N, 0.0);
    out.obsXYZ_TElocal.assign(nObs * 3, 0.0);

    // Derived freestream velocity (same convention as run_forward.cpp).
    const Real Uinf = Real(Re * nu / chord);
    (void)Ma;  // Ma accepted for API symmetry; the Amiet Mach follows calc_OASPL's
               // Uinf/340 call convention below (so the _vec path reproduces the
               // fixed-size forward path bit-for-bit on a matching grid).

    // ── 1. omega = 2*pi*f ─────────────────────────────────────────────────────
    std::vector<Real> omega(N);
    for (std::size_t i = 0; i < N; ++i)
        omega[i] = Real(2.0 * M_PI) * Real(freqs_Hz[i]);

    // ── per-surface skip flags ────────────────────────────────────────────────
    auto all_zero = [](const std::vector<double>& v) {
        return std::all_of(v.begin(), v.end(), [](double x){ return x == 0.0; });
    };
    bool upper_skip, lower_skip;
    if (has_custom) {
        // Per-surface skip: an all-zero column contributes nothing.
        upper_skip = custom_upper.empty() || all_zero(custom_upper);
        lower_skip = custom_lower.empty() || all_zero(custom_lower);
    } else {
        // tau_max <= 0 (fully-laminar TE) -> skip that surface (matches sound.hpp).
        upper_skip = !(tauMax[0] > 0.0);
        lower_skip = !(tauMax[1] > 0.0);
    }

    // ── 2. WPS stage ──────────────────────────────────────────────────────────
    std::vector<Real> WPS_upper(N, Real(0.0));
    std::vector<Real> WPS_lower(N, Real(0.0));

    if (has_custom) {
        if (!upper_skip)
            for (std::size_t i = 0; i < N; ++i) WPS_upper[i] = Real(custom_upper[i]);
        if (!lower_skip)
            for (std::size_t i = 0; i < N; ++i) WPS_lower[i] = Real(custom_lower[i]);
    } else {
        if (!upper_skip) {
            Real tauW_u = std::abs(Real(tauWall[0]));   // abs first, as existing code does
            calc_WPS_vec<Real>(model,
                               Real(theta[0]), Real(deltaStar[0]), Real(delta99[0]),
                               tauW_u, Real(tauMax[0]), Real(Ue[0]), Real(dpdx[0]),
                               omega, Real(nu), Uinf, Real(span), Real(rho),
                               1 /*isSuction: upper*/, WPS_upper);
        }
        if (!lower_skip) {
            Real tauW_l = std::abs(Real(tauWall[1]));
            calc_WPS_vec<Real>(model,
                               Real(theta[1]), Real(deltaStar[1]), Real(delta99[1]),
                               tauW_l, Real(tauMax[1]), Real(Ue[1]), Real(dpdx[1]),
                               omega, Real(nu), Uinf, Real(span), Real(rho),
                               0 /*isSuction: lower*/, WPS_lower);
        }
    }

    for (std::size_t i = 0; i < N; ++i) {
        out.WPS_upper[i] = WPS_upper[i].getValue();
        out.WPS_lower[i] = WPS_lower[i].getValue();
    }

    // ── 3. Amiet stage (per observer) ─────────────────────────────────────────
    // Edge velocities from the BL inputs; a skipped surface contributes nothing
    // regardless (its WPS column is zeros), so pass Uinf as the fallback Ue —
    // mirroring calc_OASPL's edgeVel = Uinf fallback.
    const Real Ue_top = upper_skip ? Uinf : Real(Ue[0]);
    const Real Ue_bot = lower_skip ? Uinf : Real(Ue[1]);

    const Real alpha_rad = Real(alphaDeg * M_PI / 180.0);
    const Real cos_a     = std::cos(alpha_rad);
    const Real sin_a     = std::sin(alpha_rad);
    const Real te_offset = Real(0.75) * Real(chord);

    // Amiet Mach follows the calc_OASPL / run_forward call convention.
    const Real M_amiet = Uinf / Real(340.0);

    for (std::size_t iObs = 0; iObs < nObs; ++iObs) {
        Real x_loc = Real(obsX[iObs]) * cos_a - Real(obsZ[iObs]) * sin_a - te_offset;
        Real y_loc = Real(obsY[iObs]);
        Real z_loc = Real(obsX[iObs]) * sin_a + Real(obsZ[iObs]) * cos_a;

        out.obsXYZ_TElocal[iObs * 3 + 0] = x_loc.getValue();
        out.obsXYZ_TElocal[iObs * 3 + 1] = y_loc.getValue();
        out.obsXYZ_TElocal[iObs * 3 + 2] = z_loc.getValue();

        std::vector<Real> ff(N, Real(0.0));
        TE_noise_outer_vec<Real>(M_amiet, Uinf, x_loc, y_loc, z_loc,
                                 Real(chord / 2.0), Real(chord), Real(span),
                                 Real(340.0), omega,
                                 Ue_bot, Ue_top,
                                 WPS_lower, WPS_upper, ff);

        for (std::size_t i = 0; i < N; ++i)
            out.FF_spectra[iObs * N + i] = ff[i].getValue();
    }

    return out;
}
