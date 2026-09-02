#pragma once

// Trailing-edge noise model — Roger & Moreau (2005). TWO paths, on purpose:
//
//  * fixed-Nsound templates (AD-critical, used by sound.hpp/calc_OASPL) — the
//    MID-SPAN model:
//        x₂ = 0, K̄₂ = 0, κ̄ = μ̄, kbar² = μ̄² > 0 → supercritical branch only.
//  * the *_vec overloads (noise_run only, never taped) — the GENERAL 3-D
//    oblique-gust model: reads x₂, forms K̄₂ = k̄x₂/S0, branches on criticality,
//    and regularises the cut. See the "General oblique-gust kernel" block near
//    the bottom of this file.
//
// The two agree exactly at x₂ = 0 (regression-gated). They diverge off mid-span
// — deliberately: the optimisation path's observers are mid-span, and restoring
// the general kernel there would change the taped tape. See CHANGELOG.
//
// Key functions:
//   Estar()                  — E*(x), R&M Eq. 3-4 (via complex erf/Faddeeva)
//   Radiation_integral1()    — I1, R&M Eq. 13 (primary TE scattering)
//   Radiation_integral2()    — I2, R&M Eq. 14 (leading-edge back-scatter)
//   Radiation_integral_total() — frequency loop over Nsound points
//   TE_noise_outer()         — far-field PSD, R&M Eq. 18
//
// CoDiPack note:
//   errFunc() uses .getValue()/.getGradient() and StatementPushHelper.
//   It requires Real to be a CoDiPack active type — will not compile
//   with Real = double. All call paths go through codi::RealReverse.

#include <cmath>
#include <complex>
#include <vector>
#include <codi.hpp>  // for codi::StatementPushHelper in errFunc
#include "Faddeeva.hh"



// divide (ar + i ai)/(br + i bi)
template<typename Real>
inline void cdiv(Real ar, Real ai, Real br, Real bi, Real &cr, Real &ci) {
    Real den = br*br + bi*bi;
    cr = (ar*br + ai*bi)/den;
    ci = (ai*br - ar*bi)/den;
}

// sqrt(ar + i ai)
template<typename Real>
inline void complex_sqrt(Real ar, Real ai,
                         Real &br, Real &bi)
{
    if (ai == 0.0) {
        if (ar >= 0.0) {
            br = std::sqrt(ar);
            bi = 0.0;
        } else {
            br = 0.0;
            bi = std::sqrt(-ar);
        }
        return;
    }

    Real r = std::sqrt(ar*ar + ai*ai);
    Real t = std::sqrt(0.5*(r + std::abs(ar)));

    if (ar >= 0.0) {
        br = t;
        bi = ai/(2.0*t);
    } else {
        br = std::abs(ai)/(2.0*t);
        bi = (ai >= 0.0) ? t : -t;
    }
}


template<typename Real>
void errFunc(Real in_r, Real in_i, Real &out_r, Real &out_i){

    double x_val = in_r.getValue();
    double y_val = in_i.getValue();
    std::complex<double> z(x_val, y_val);
    std::complex<double> w = Faddeeva::erf(z);

    std::complex<double> dw_dz = (2.0 / std::sqrt(M_PI)) * std::exp(-z * z);

    double du_dx = dw_dz.real();
    double du_dy = -dw_dz.imag();
    double dv_dx = dw_dz.imag();
    double dv_dy = dw_dz.real();

    codi::StatementPushHelper<Real> ph;
    ph.startPushStatement();
    ph.pushArgument(in_r, du_dx);
    ph.pushArgument(in_i, du_dy);
    ph.endPushStatement(out_r, w.real());

    codi::StatementPushHelper<Real> phIm;
    phIm.startPushStatement();
    phIm.pushArgument(in_r, dv_dx);
    phIm.pushArgument(in_i, dv_dy);
    phIm.endPushStatement(out_i, w.imag());
}


// erfcx(z) = e^{z²}·erfc(z) — the SCALED complex error function, via Faddeeva.
// Same CoDiPack external-function pattern as errFunc above (and the same
// restriction: Real must be a CoDiPack active type).
//   d/dz erfcx(z) = 2z·erfcx(z) − 2/√π
// Needed only by the subcritical branch of the general kernel, where the exact
// product e^{−2iA′₁}·erf(ζ) is O(1) but its factors are e^{∓2κ̄′}: forming them
// separately overflows to inf×0 = NaN for κ̄′ ≳ 350. Since ζ² = −2iA′₁ exactly,
//   e^{−2iA′₁}·erf(ζ) = e^{ζ²}(1 − erfc(ζ)) = e^{ζ²} − erfcx(ζ),
// and both terms on the right stay bounded.
template<typename Real>
void erfcxFunc(Real in_r, Real in_i, Real &out_r, Real &out_i){

    std::complex<double> z(in_r.getValue(), in_i.getValue());
    std::complex<double> w = Faddeeva::erfcx(z);

    std::complex<double> dw_dz = 2.0 * z * w - 2.0 / std::sqrt(M_PI);

    double du_dx =  dw_dz.real();
    double du_dy = -dw_dz.imag();
    double dv_dx =  dw_dz.imag();
    double dv_dy =  dw_dz.real();

    codi::StatementPushHelper<Real> ph;
    ph.startPushStatement();
    ph.pushArgument(in_r, du_dx);
    ph.pushArgument(in_i, du_dy);
    ph.endPushStatement(out_r, w.real());

    codi::StatementPushHelper<Real> phIm;
    phIm.startPushStatement();
    phIm.pushArgument(in_r, dv_dx);
    phIm.pushArgument(in_i, dv_dy);
    phIm.endPushStatement(out_i, w.imag());
}


/////////////////////////////////////////// Amiet model — Roger & Moreau (2005) ///////////////////

// Estar(xr + i·xi) — Roger & Moreau (2005) Eq. 3-4.
// Computes E*(x) = erf[(1+i)·sqrt((xr+i·xi)/2)] / (1+i)
// For real x (xi=0): mid-span E* used in radiation integrals.
// The no-star E(x) = conj(E*(x)) for real x — obtained by negating
// the imaginary part of Estar at the call site.
template<typename Real>
inline void Estar(Real xr, Real xi,
                  Real &Er, Real &Ei)
{
    Real sr, si;
    complex_sqrt<Real>(0.5*xr, 0.5*xi, sr, si);
    Real mr = sr - si;   // Re[(1+i)·sqrt(x/2)]
    Real mi = sr + si;   // Im[(1+i)·sqrt(x/2)]

    Real erf_r, erf_i;
    errFunc<Real>(mr, mi, erf_r, erf_i);

    // erf(z) / (1+i) = ((a+b) + i(b-a)) / 2
    Er = 0.5*(erf_r + erf_i);
    Ei = 0.5*(erf_i - erf_r);
}

// sinc(z) = sin(z)/z, limit 1 as z→0; used for G_a/G_b in Radiation_integral2
template<typename Real>
inline Real sinc_safe(Real z) {
    return (std::abs(z) < Real(1e-10)) ? Real(1.0) : std::sin(z) / z;
}

// Primary TE scattering term, Roger & Moreau (2005) Eq.13. Supercritical only
// (K2_bar=0 ⟹ kappa_bar=mu_bar always). B = K̄₁ + (1+M)μ̄ (at mid-span), C from Eq.12.
template<typename Real>
inline void Radiation_integral1(Real B, Real C,
                                Real &f1r, Real &f1i)
{
    Real a_r, a_i; Estar<Real>(2.0*(B-C), 0.0, a_r, a_i);
    Real b_r, b_i; Estar<Real>(2.0*B,     0.0, b_r, b_i);

    Real cos2C = std::cos(2.0*C), sin2C = std::sin(2.0*C);
    // prefactor = -e^{2iC}/(iC) = i·e^{2iC}/C
    // Re = -sin(2C)/C,  Im = +cos(2C)/C
    Real pref_r = -sin2C / C;
    Real pref_i =  cos2C / C;

    Real onepI_r = 1.0, onepI_i = 1.0;  // (1+i)
    Real e_2C_r = std::cos(-2.0*C), e_2C_i = std::sin(-2.0*C);
    Real s = std::sqrt(2.0*B);

    Real sc_r, sc_i;
    complex_sqrt<Real>(2.0*(B-C), 0.0, sc_r, sc_i);

    Real tmp_r = onepI_r*e_2C_r - onepI_i*e_2C_i;
    Real tmp_i = onepI_r*e_2C_i + onepI_i*e_2C_r;
    tmp_r *= s; tmp_i *= s;

    Real a_div_r, a_div_i;
    cdiv<Real>(a_r, a_i, sc_r, sc_i, a_div_r, a_div_i);

    // Combine s and a_div: sqrt(2B) * E*(2(B-C))/sqrt(2(B-C))
    //                    = sqrt(B/(B-C)) * E*(2(B-C))  per R&M Eq.13
    Real t1r = tmp_r*a_div_r - tmp_i*a_div_i;
    Real t1i = tmp_r*a_div_i + tmp_i*a_div_r;

    Real t2r = -(onepI_r*b_r - onepI_i*b_i);
    Real t2i = -(onepI_r*b_i + onepI_i*b_r);

    Real br_r = t1r + t2r + 1.0;
    Real br_i = t1i + t2i;

    f1r = pref_r*br_r - pref_i*br_i;
    f1i = pref_r*br_i + pref_i*br_r;
}


// Roger & Moreau (2005) Eq. 14 — leading-edge back-scattering correction I2.
// Supercritical mid-span: kappa_bar = mu_bar, k_min_bar = mu_bar.
// H = prefactor (R&M Eq.14 notation), epsilon = (1 + 1/(4*mu_bar))^(-1/2) (Eq.9).
// G = sum of sub-terms G_a..G_e; Estar(x) for E*, conj(Estar(x)) for E(x).
// NOTE: G_e coefficient uses sqrt(0.5*k/D), corrected from sqrt(k/D).
// TODO: replace denominator clipping in G_c, G_d, G_e with limiting forms.
template<typename Real>
void Radiation_integral2(
    Real B, Real K_bar, Real k_min_bar, Real mu_bar, Real S0,
    Real K_1_bar, Real alpha, Real x, Real M,
    Real &f2r, Real &f2i)
{
    Real error = std::pow(1.0 + 1.0/(4.0*mu_bar), -0.5);   // Eq.9
    // D = mu_bar*(1 - x/S0); k_min_bar = mu_bar at mid-span (K_2_bar=0)
    Real D = mu_bar * (1.0 - x / S0);

    // Cache Estar(4·k_min_bar) — reused in Ẽ block, G_c, and G_d.
    Real Fr_4k, Fi_4k;
    Estar<Real>(4.0*k_min_bar, 0.0, Fr_4k, Fi_4k);

    // Ẽ = exp(4i·k_min_bar)·(1 - (1+i)·E*(4·k_min_bar))
    Real Fr, Fi;
    Fr = Fr_4k; Fi = Fi_4k;
    Real t1r = Fr - Fi;   // (1+i)·E*: real part
    Real t1i = Fi + Fr;   // (1+i)·E*: imag part
    Real oneMinus_r = 1.0 - t1r;
    Real oneMinus_i =     - t1i;
    Real e4r = std::cos(4.0*k_min_bar);
    Real e4i = std::sin(4.0*k_min_bar);
    Real Er = e4r*oneMinus_r - e4i*oneMinus_i;
    Real Ei = e4r*oneMinus_i + e4i*oneMinus_r;
    // imaginary-part correction: Ẽ = Re(E) + i·ε·Im(E)
    Real Efr = Er;
    Real Efi = error*Ei;

    // --- G_a: (1+ε)·e^{i(2k+D)}·sinc(D-2k)
    Real phase = 2.0*k_min_bar + D;
    Real epr = std::cos(phase);
    Real epi = std::sin(phase);
    Real sinc_a = sinc_safe<Real>(D - 2.0*k_min_bar);
    Real G_ar = (1.0+error)*epr*sinc_a;
    Real G_ai = (1.0+error)*epi*sinc_a;

    // --- G_b: (1-ε)·e^{i(-2k+D)}·sinc(D+2k)
    phase = -2.0*k_min_bar + D;
    epr = std::cos(phase);
    epi = std::sin(phase);
    Real sinc_b = sinc_safe<Real>(D + 2.0*k_min_bar);
    Real G_br = (1.0-error)*epr*sinc_b;
    Real G_bi = (1.0-error)*epi*sinc_b;

    // --- G_c: [(1+ε)(1-i)] / [2(D-2k)] · e^{4ik}·E*(4k)
    Real denC_val = D - 2.0*k_min_bar;
    // TODO: replace denominator clipping with analytical limiting form near
    // D +/- 2*k_min_bar = 0. Current clipping prevents crashes but may
    // introduce gradient discontinuities for CoDiPack AD.
    Real denC = 2.0 * ((std::abs(denC_val) < Real(1e-10)) ? Real(1e-10) : denC_val);
    Real m1r = 1.0, m1i = -1.0;  // (1-i)
    Real coeffr = (1.0+error)*m1r / denC;
    Real coeffi = (1.0+error)*m1i / denC;
    epr = std::cos(4.0*k_min_bar);
    epi = std::sin(4.0*k_min_bar);
    Fr = Fr_4k; Fi = Fi_4k;
    Real tmp_r = epr*Fr - epi*Fi;
    Real tmp_i = epr*Fi + epi*Fr;
    Real G_cr = coeffr*tmp_r - coeffi*tmp_i;
    Real G_ci = coeffr*tmp_i + coeffi*tmp_r;

    // --- G_d: [(1-ε)(1+i)] / [2(D+2k)] · e^{-4ik}·E(4k),  subtracted in sum
    Real denD_val = D + 2.0*k_min_bar;
    // TODO: replace denominator clipping with analytical limiting form near
    // D +/- 2*k_min_bar = 0. Current clipping prevents crashes but may
    // introduce gradient discontinuities for CoDiPack AD.
    Real denD = 2.0 * ((std::abs(denD_val) < Real(1e-10)) ? Real(1e-10) : denD_val);
    Real p1r = 1.0, p1i = 1.0;  // (1+i)
    Real coeffDr = (1.0-error)*p1r / denD;
    Real coeffDi = (1.0-error)*p1i / denD;
    epr = std::cos(-4.0*k_min_bar);
    epi = std::sin(-4.0*k_min_bar);
    // E(4k) = conj(E*(4k)): negate imaginary part of cached E*(4k)
    Fr = Fr_4k; Fi = -Fi_4k;
    tmp_r = epr*Fr - epi*Fi;
    tmp_i = epr*Fi + epi*Fr;
    Real G_dr = coeffDr*tmp_r - coeffDi*tmp_i;
    Real G_di = coeffDr*tmp_i + coeffDi*tmp_r;

    // --- G_e: [e^{2iD}/2]·sqrt(2k/D)·E*(2D)·bracket; guard D=0
    Real G_er, G_ei;
    if (std::abs(D) < Real(1e-10)) {
        G_er = 0.0; G_ei = 0.0;
    } else {
        Real e2r = std::cos(2.0*D);
        Real e2i = std::sin(2.0*D);
        // R&M Eq.14: e^{2iD}/2 * sqrt(2k/D) = e^{2iD} * sqrt(k/(2D))
        // so the sqrt argument is 0.5*k_min_bar/D, not k_min_bar/D.
        Real sqrtfactor_r, sqrtfactor_i;
        complex_sqrt<Real>(0.5*k_min_bar / D, Real(0.0), sqrtfactor_r, sqrtfactor_i);

        Estar<Real>(2.0*D, 0.0, Fr, Fi);
        tmp_r = e2r*Fr - e2i*Fi;
        tmp_i = e2r*Fi + e2i*Fr;
        Real new_r = tmp_r*sqrtfactor_r - tmp_i*sqrtfactor_i;
        Real new_i = tmp_r*sqrtfactor_i + tmp_i*sqrtfactor_r;
        tmp_r = new_r; tmp_i = new_i;

        // bracket = (1+i)·(1-ε)/(D+2k) - (1-i)·(1+ε)/(D-2k)
        Real da = D + 2.0*k_min_bar;
        Real db = D - 2.0*k_min_bar;
        // TODO: replace clipping with limiting forms near D +/- 2k = 0.
        Real da_s = (std::abs(da) < Real(1e-10)) ? Real(1e-10) : da;
        Real db_s = (std::abs(db) < Real(1e-10)) ? Real(1e-10) : db;
        Real term1r = (1.0-error) / da_s;
        Real term1i =  term1r;    // (1+i)·coeff: r=i=coeff
        Real term2r = (1.0+error) / db_s;
        Real term2i = -term2r;    // (1-i)·coeff: r=coeff, i=-coeff
        Real Br_r = term1r - term2r;
        Real Br_i = term1i - term2i;
        G_er = tmp_r*Br_r - tmp_i*Br_i;
        G_ei = tmp_r*Br_i + tmp_i*Br_r;
    }

    // G_d subtracted: spec has minus sign for the E(4k) term
    Real G_r = G_ar + G_br + G_cr - G_dr + G_er;
    Real G_i = G_ai + G_bi + G_ci - G_di + G_ei;

    // --- H = (1+i)·e^{-4ik} / [2√π·(α-1)·K̄·sqrt(B)] · (1-Y²)
    Real Theta = std::sqrt((K_1_bar + (1.0+M)*mu_bar) / (K_bar + (1.0+M)*mu_bar));
    Real Hcoeff = (1.0 - Theta*Theta) / (2.0*std::sqrt(M_PI)*(alpha-1.0)*K_bar*std::sqrt(B));
    epr = std::cos(-4.0*k_min_bar);
    epi = std::sin(-4.0*k_min_bar);
    Real m1pr = 1.0, m1pi = 1.0;  // (1+i)
    Real Hr = m1pr*epr - m1pi*epi;
    Real Hi = m1pr*epi + m1pi*epr;
    Hr *= Hcoeff; Hi *= Hcoeff;

    // f2 = H·(Ẽ - e^{2iD} + i·(D + K̄ + (M-1)·μ̄)·G)
    epr = std::cos(2.0*D);
    epi = std::sin(2.0*D);
    Real part1r = Efr - epr;
    Real part1i = Efi - epi;

    // i·coeffI·G = -coeffI·G_i + i·coeffI·G_r
    Real coeffI = D + K_bar + (M - 1.0)*mu_bar;
    Real part2r = -coeffI*G_i;
    Real part2i =  coeffI*G_r;

    Real totalr = part1r + part2r;
    Real totali = part1i + part2i;

    f2r = Hr*totalr - Hi*totali;
    f2i = Hr*totali + Hi*totalr;
}

// Loops over Nsound frequencies, assembles I_abs2[Nsound] for far-field PSD.
// R&M (2005) Section 3.1 — mid-span, K̄₂=0 always ⟹ always supercritical.
template<typename Real>
void Radiation_integral_total(
    const Real *C,
    const Real *K_bar,
    const Real *mu_bar,
    Real S0,
    const Real *K_1_bar,
    Real alpha,
    Real M,
    Real x,
    Real *I_abs2
)
{
    for (int i = 0; i < Nsound; ++i)
    {
        // K̄₂ = 0 at mid-span → kbar2 = μ̄², k_min_bar = μ̄, B = K̄₁ + (1+M)·μ̄
        Real k_min_bar = mu_bar[i];
        Real B = K_1_bar[i] + M*mu_bar[i] + k_min_bar;

        Real fr1, fi1;
        Radiation_integral1<Real>(B, C[i], fr1, fi1);

        Real fr2, fi2;
        Radiation_integral2<Real>(
            B,
            K_bar[i],
            k_min_bar,
            mu_bar[i],
            S0,
            K_1_bar[i],
            alpha,
            x,
            M,
            fr2, fi2);

        Real Ireal = fr1 + fr2;
        Real Iimag = fi1 + fi2;

        I_abs2[i] = Ireal*Ireal + Iimag*Iimag;
    }
}


// Roger & Moreau (2005) Eq. 18 — far-field PSD S_pp(omega).
// Mid-span observer (y=0): S0 = sqrt(x^2 + beta^2*z^2).
// Integrates upper and lower surface contributions.
// Inputs:
//   M, U    — freestream Mach and velocity
//   x, z    — observer coordinates (y ignored, assumed 0 for mid-span)
//   b       — semi-chord (half-chord = c/2), used for K_bar = omega/U * b
//   c       — full chord (used only for b_half = b, kept for clarity)
//   span    — span length
//   c0      — speed of sound
//   omega   — angular frequency array [Nsound]
//   Ue_bot/top — trailing-edge edge velocity, lower/upper surface
//   WPS_lower/upper — one-sided wall-pressure spectrum Phi_pp [Nsound]
//   farfieldSpectra — output: S_pp(omega) [Nsound]
template<typename Real>
void TE_noise_outer(
    Real M, Real U, Real x, Real y, Real z,
    Real b, Real c, Real span,

    Real c0,

    const Real omega[Nsound],
    const Real Ue_bot, const Real Ue_top,
    Real (&WPS_lower)[Nsound], Real (&WPS_upper)[Nsound],

    Real (&farfieldSpectra)[Nsound]
)
{
    Real Ue[2] = {Ue_top, Ue_bot};

    // beta and S0 depend only on M, x, z — hoist outside the surf loop.
    // Mid-span implementation: y is ignored and assumed to be zero.
    // Including y in S0 while K2_bar=0 is physically inconsistent
    // (R&M Section 3, mid-span observer: x2=0 throughout).
    Real beta = std::sqrt(1.0 - M*M);
    Real S0   = std::sqrt(x*x + beta*beta*z*z);

    for (int surf = 0; surf < 2; ++surf) {

    Real U_c  = 0.7 * Ue[surf];
    Real alpha = U / U_c;

    // Roger & Moreau (2005) Section 2 — non-dimensional wavenumbers.
    // K̄₂ = 0 (mid-span, y=0); kbar² = μ̄² always supercritical.
    Real C[Nsound];
    Real K_bar[Nsound];
    Real mu_bar[Nsound];
    Real K_1_bar[Nsound];

    for (int i = 0; i < Nsound; ++i)
    {
        Real K = omega[i] / U;

        K_bar[i]   = K * b;
        mu_bar[i]  = K_bar[i] * M / (beta*beta);
        K_1_bar[i] = alpha * K_bar[i];
        C[i]       = K_1_bar[i] - mu_bar[i] * (x/S0 - M);   // Eq. 12
    }

    Real I_abs2[Nsound];
    Radiation_integral_total(C, K_bar, mu_bar, S0, K_1_bar, alpha, M, x, I_abs2);

    // Roger & Moreau (2005) Eq. 19 — Corcos spanwise correlation length (K̄₂=0).
    Real b_c = 1.47;
    Real l_y[Nsound];
    for (int i = 0; i < Nsound; ++i)
        l_y[i] = (b_c * U_c) / omega[i];

    // Roger & Moreau (2005) Eq. 18 — far-field PSD S_pp(ω).
    // S_pp = (ωb·z / 2πc₀S₀²)² · 2·span · |I|² · Φ_pp · l_y,
    // b = semi-chord (half-chord), passed directly as parameter
    Real b_half = b;  // b is already the semi-chord (c/2) passed by caller
    for (int i = 0; i < Nsound; ++i) {
        Real _t1   = (omega[i]*b_half*z) / (2.0*M_PI*c0*S0*S0);
        Real term1 = _t1 * _t1;

        if (surf == 0) {
            farfieldSpectra[i]  = term1 * 2.0*span * I_abs2[i] * WPS_upper[i] * l_y[i];
        } else {
            farfieldSpectra[i] += term1 * 2.0*span * I_abs2[i] * WPS_lower[i] * l_y[i];
        }
    }
    }

}


/////////////////////////// Length-generic (_vec) overloads ///////////////////////////
//
// The fixed [Nsound] arrays become std::vector<Real> of length omega.size().
// Used ONLY by noise_run (acoustics-only forward path, never AD'd); the
// fixed-Nsound versions above remain on the AD-critical path untouched.
//
// PHYSICS DIVERGENCE (phase 2, July 2026). Unlike the fixed-Nsound path above,
// which is the Roger & Moreau (2005) mid-span (K̄₂ = 0) model, the _vec path
// implements the GENERAL three-dimensional oblique-gust formulation: it reads
// the spanwise observer coordinate x₂, selects the gust of Eq. 18, and takes
// the supercritical or subcritical branch as the criticality parameter demands.
// The two agree exactly at x₂ = 0. See the header block below.


/////////////////// General oblique-gust (K̄₂ ≠ 0) kernel — _vec only ///////////////////
//
// Restores the general 3-D gust formulation of Roger & Moreau (2005) that the
// May 2026 audit removed as unused. Confined to the _vec overloads: the
// fixed-Nsound templates above are the K̄₂ = 0 instance of the same result and
// stay byte-identical on the AD-critical path.
//
// Notation (the paper's; overbar = non-dimensionalised by the semi-chord b):
//   S0² = x₁² + β²(x₂² + x₃²)        general convected source-observer distance
//   K̄₂  = k̄·x₂/S0                   the gust Eq. 18's delta selection picks out
//   ξ   = K̄₂/(βμ̄)                   criticality parameter
//   κ̄   = μ̄·sqrt(1 − ξ²)            supercritical scattered wavenumber (ξ < 1)
//   κ̄′  = μ̄·sqrt(ξ² − 1)            subcritical decay rate            (ξ > 1)
// equivalently κ̄² = μ̄² − K̄₂²/β².
//
// IDENTITY:  ξ = β|x₂|/S0.
// Proof: k̄ = KbM = K̄M and μ̄ = K̄M/β², so μ̄ = k̄/β². Hence
//     ξ = K̄₂/(βμ̄) = (k̄ x₂/S0)/(β·k̄/β²) = β x₂/S0.
// Since S0² = x₁² + β²x₂² + β²x₃² ≥ β²x₂², we get ξ ≤ 1 ALWAYS, with equality
// only at x₁ = x₃ = 0. Two consequences:
//   * the gust selected by Eq. 18 is always supercritical in production
//     (noise_run); it approaches the cut κ̄ → 0 only as the observer nears the
//     blade's spanwise axis — exactly the geometry the mid-span kernel got
//     wrong (see CHANGELOG, "Rotor TE noise post-processing", phase-1
//     limitation);
//   * the subcritical forms are still REQUIRED — as the far-side anchor of the
//     near-cutoff bridge below, and for validation sweeps that drive K̄₂ as a
//     free parameter across the cut (as the paper's Figs. 9-11 do).
//
// Only K̄₂² enters |I|², l_y and κ̄, so the sign of x₂ never matters; the signed
// value is nonetheless carried through.

// Near-cutoff regularisation half-width, in κ̄. R&M §4.1 attribute the deep,
// narrow cuts at κ̄ = 0 to non-convergence of the two-step Schwarzschild
// iteration there, and state the back-scattering approximation is "expected to
// be accurate enough for κ̄ > 0.125". That sentence is this constant's sole
// provenance.
constexpr double AMIET_KAPPA_REG = 0.125;


// R&M Eq. 19 — Corcos spanwise correlation length, general spanwise wavenumber:
//     l_y(K₂,ω) = (ω/(b_c·U_c)) / (K₂² + (ω/(b_c·U_c))²)
// Evaluated in the algebraically identical form
//     l_y = l_c / (1 + (K₂·l_c)²),   l_c = b_c·U_c/ω,
// which (a) reduces at K₂ = 0 to exactly l_c — the mid-span expression divided
// by 1.0, hence bit-identical, not merely equal to round-off — and (b) never
// forms ω/(b_c·U_c) explicitly, so it cannot overflow at large ω.
// Monotone decreasing in |K₂| and strictly positive for U_c, ω > 0.
template<typename Real>
inline Real corcos_l_y(Real K_2, Real omega, Real b_c, Real U_c)
{
    Real l_c = (b_c * U_c) / omega;
    Real q   = K_2 * l_c;
    return l_c / (1.0 + q*q);
}


// R&M Eq. 13 (supercritical, B real) AND Eq. 15 (subcritical, B complex).
//
// Eq. 15 is the EXACT analytic continuation κ̄ → −iκ̄′ of Eq. 13:
//     B   = K̄₁ + Mμ̄ + κ̄     ->  A′₁ = K̄₁ + Mμ̄ − iκ̄′
//     B−C = κ̄ + μ̄·x₁/S0      ->  Z   = μ̄·x₁/S0 − iκ̄′
// (the second follows from the first because C = K̄₁ − μ̄(x₁/S0 − M) carries no
// K̄₂ and is real on both branches). So one complex-B implementation serves
// both, and I₁ is continuous across the cut by construction. Writing Eq. 15
// with the paper's F* and applying F*(sqrt(i·w)) = (1+i)·E*(w) reproduces this
// expression term for term.
//
//   I₁ = i·e^{2iC}/C · { (1+i)·e^{−2iC}·sqrt(B/(B−C))·E*(2(B−C))
//                        − (1+i)·E*(2B) + 1 }
//
// At Bi = 0 this reduces expression-for-expression (and bit-for-bit) to
// Radiation_integral1 above, which is already general in (B, C).
template<typename Real>
void Radiation_integral1_cplx(Real Br, Real Bi, Real C,
                              Real &f1r, Real &f1i)
{
    Real BmCr = Br - C, BmCi = Bi;

    Real a_r, a_i; Estar<Real>(2.0*BmCr, 2.0*BmCi, a_r, a_i);
    Real b_r, b_i; Estar<Real>(2.0*Br,   2.0*Bi,   b_r, b_i);

    Real cos2C = std::cos(2.0*C), sin2C = std::sin(2.0*C);
    // prefactor = -e^{2iC}/(iC) = i·e^{2iC}/C
    Real pref_r = -sin2C / C;
    Real pref_i =  cos2C / C;

    Real onepI_r = 1.0, onepI_i = 1.0;  // (1+i)
    Real e_2C_r = std::cos(-2.0*C), e_2C_i = std::sin(-2.0*C);

    Real s_r, s_i;   complex_sqrt<Real>(2.0*Br,   2.0*Bi,   s_r, s_i);
    Real sc_r, sc_i; complex_sqrt<Real>(2.0*BmCr, 2.0*BmCi, sc_r, sc_i);

    // B−C = κ̄ + μ̄x₁/S0 can approach zero near the cut when x₁ < 0. The ratio
    // E*(2z)/sqrt(2z) has the finite limit sqrt(2/π) as z→0, but the quotient
    // itself is 0/0; clip the denominator, matching the guard style used in
    // Radiation_integral2 above.
    Real sc_mod = std::sqrt(sc_r*sc_r + sc_i*sc_i);
    if (sc_mod < Real(1e-10)) { sc_r = Real(1e-10); sc_i = Real(0.0); }

    Real tmp_r = onepI_r*e_2C_r - onepI_i*e_2C_i;
    Real tmp_i = onepI_r*e_2C_i + onepI_i*e_2C_r;
    // multiply by sqrt(2B)
    Real ts_r = tmp_r*s_r - tmp_i*s_i;
    Real ts_i = tmp_r*s_i + tmp_i*s_r;

    Real a_div_r, a_div_i;
    cdiv<Real>(a_r, a_i, sc_r, sc_i, a_div_r, a_div_i);

    Real t1r = ts_r*a_div_r - ts_i*a_div_i;
    Real t1i = ts_r*a_div_i + ts_i*a_div_r;

    Real t2r = -(onepI_r*b_r - onepI_i*b_i);
    Real t2i = -(onepI_r*b_i + onepI_i*b_r);

    Real br_r = t1r + t2r + 1.0;
    Real br_i = t1i + t2i;

    f1r = pref_r*br_r - pref_i*br_i;
    f1i = pref_r*br_i + pref_i*br_r;
}


// R&M Eq. 14 — leading-edge back-scattering correction, GENERAL κ̄.
//
// This is Radiation_integral2 above with the substitution map of the phase-2
// implementation applied: every additive μ̄ acting as the equivalent 2-D
// frequency parameter becomes κ̄, while the μ̄ arising from the change of
// variables (the Mμ̄ convection terms and the μ̄·x₁/S0 phase) stays μ̄. Four
// expressions change; each reduces identically at κ̄ = μ̄:
//     ε      : argument μ̄ -> κ̄        (functional form preserved, see below)
//     D      : μ̄(1 − x₁/S0)  ->  κ̄ − μ̄·x₁/S0
//     Y²     : (K̄₁+(1+M)μ̄)/(K̄+(1+M)μ̄)  ->  (K̄₁+Mμ̄+κ̄)/(K̄+Mμ̄+κ̄) = B/A
//     coeffI : D + K̄ + (M−1)μ̄  ->  D + K̄ + Mμ̄ − κ̄
// Everything else already carried kappa_bar as a separate parameter.
//
// ε FORM: kept exactly as the validated mid-span code has it, (1+1/(4κ̄))^(−1/2),
// with only the argument swapped. The project summary document renders Eq. 9 as
// 1 + (4μ̄)^(−1/2) instead; that reading is NOT adopted here and is flagged in
// CHANGELOG for a check against the typeset paper.
template<typename Real>
void Radiation_integral2_general(
    Real B, Real K_bar, Real kappa_bar, Real mu_bar, Real S0,
    Real K_1_bar, Real alpha, Real x, Real M,
    Real &f2r, Real &f2i)
{
    Real error = std::pow(1.0 + 1.0/(4.0*kappa_bar), -0.5);   // Eq.9, arg κ̄
    Real D = kappa_bar - mu_bar * x / S0;
    Real k_min_bar = kappa_bar;

    Real Fr_4k, Fi_4k;
    Estar<Real>(4.0*k_min_bar, 0.0, Fr_4k, Fi_4k);

    // Ẽ = exp(4i·κ̄)·(1 - (1+i)·E*(4κ̄)),  then Ẽ = Re(Ẽ) + i·ε·Im(Ẽ)
    Real Fr = Fr_4k, Fi = Fi_4k;
    Real t1r = Fr - Fi;
    Real t1i = Fi + Fr;
    Real oneMinus_r = 1.0 - t1r;
    Real oneMinus_i =     - t1i;
    Real e4r = std::cos(4.0*k_min_bar);
    Real e4i = std::sin(4.0*k_min_bar);
    Real Er = e4r*oneMinus_r - e4i*oneMinus_i;
    Real Ei = e4r*oneMinus_i + e4i*oneMinus_r;
    Real Efr = Er;
    Real Efi = error*Ei;

    // --- G_a: (1+ε)·e^{i(2κ̄+D)}·sinc(D-2κ̄)
    Real phase = 2.0*k_min_bar + D;
    Real epr = std::cos(phase);
    Real epi = std::sin(phase);
    Real sinc_a = sinc_safe<Real>(D - 2.0*k_min_bar);
    Real G_ar = (1.0+error)*epr*sinc_a;
    Real G_ai = (1.0+error)*epi*sinc_a;

    // --- G_b: (1-ε)·e^{i(-2κ̄+D)}·sinc(D+2κ̄)
    phase = -2.0*k_min_bar + D;
    epr = std::cos(phase);
    epi = std::sin(phase);
    Real sinc_b = sinc_safe<Real>(D + 2.0*k_min_bar);
    Real G_br = (1.0-error)*epr*sinc_b;
    Real G_bi = (1.0-error)*epi*sinc_b;

    // --- G_c: [(1+ε)(1-i)] / [2(D-2κ̄)] · e^{4iκ̄}·E*(4κ̄)
    Real denC_val = D - 2.0*k_min_bar;
    Real denC = 2.0 * ((std::abs(denC_val) < Real(1e-10)) ? Real(1e-10) : denC_val);
    Real m1r = 1.0, m1i = -1.0;
    Real coeffr = (1.0+error)*m1r / denC;
    Real coeffi = (1.0+error)*m1i / denC;
    epr = std::cos(4.0*k_min_bar);
    epi = std::sin(4.0*k_min_bar);
    Fr = Fr_4k; Fi = Fi_4k;
    Real tmp_r = epr*Fr - epi*Fi;
    Real tmp_i = epr*Fi + epi*Fr;
    Real G_cr = coeffr*tmp_r - coeffi*tmp_i;
    Real G_ci = coeffr*tmp_i + coeffi*tmp_r;

    // --- G_d: [(1-ε)(1+i)] / [2(D+2κ̄)] · e^{-4iκ̄}·E(4κ̄),  subtracted in sum
    Real denD_val = D + 2.0*k_min_bar;
    Real denD = 2.0 * ((std::abs(denD_val) < Real(1e-10)) ? Real(1e-10) : denD_val);
    Real p1r = 1.0, p1i = 1.0;
    Real coeffDr = (1.0-error)*p1r / denD;
    Real coeffDi = (1.0-error)*p1i / denD;
    epr = std::cos(-4.0*k_min_bar);
    epi = std::sin(-4.0*k_min_bar);
    Fr = Fr_4k; Fi = -Fi_4k;          // E(4κ̄) = conj(E*(4κ̄))
    tmp_r = epr*Fr - epi*Fi;
    tmp_i = epr*Fi + epi*Fr;
    Real G_dr = coeffDr*tmp_r - coeffDi*tmp_i;
    Real G_di = coeffDr*tmp_i + coeffDi*tmp_r;

    // --- G_e: [e^{2iD}/2]·sqrt(2κ̄/D)·E*(2D)·bracket; guard D=0
    Real G_er, G_ei;
    if (std::abs(D) < Real(1e-10)) {
        G_er = 0.0; G_ei = 0.0;
    } else {
        Real e2r = std::cos(2.0*D);
        Real e2i = std::sin(2.0*D);
        Real sqrtfactor_r, sqrtfactor_i;
        complex_sqrt<Real>(0.5*k_min_bar / D, Real(0.0), sqrtfactor_r, sqrtfactor_i);

        Estar<Real>(2.0*D, 0.0, Fr, Fi);
        tmp_r = e2r*Fr - e2i*Fi;
        tmp_i = e2r*Fi + e2i*Fr;
        Real new_r = tmp_r*sqrtfactor_r - tmp_i*sqrtfactor_i;
        Real new_i = tmp_r*sqrtfactor_i + tmp_i*sqrtfactor_r;
        tmp_r = new_r; tmp_i = new_i;

        Real da = D + 2.0*k_min_bar;
        Real db = D - 2.0*k_min_bar;
        Real da_s = (std::abs(da) < Real(1e-10)) ? Real(1e-10) : da;
        Real db_s = (std::abs(db) < Real(1e-10)) ? Real(1e-10) : db;
        Real term1r = (1.0-error) / da_s;
        Real term1i =  term1r;
        Real term2r = (1.0+error) / db_s;
        Real term2i = -term2r;
        Real Br_r = term1r - term2r;
        Real Br_i = term1i - term2i;
        G_er = tmp_r*Br_r - tmp_i*Br_i;
        G_ei = tmp_r*Br_i + tmp_i*Br_r;
    }

    Real G_r = G_ar + G_br + G_cr - G_dr + G_er;
    Real G_i = G_ai + G_bi + G_ci - G_di + G_ei;

    // --- H = (1+i)·e^{-4iκ̄} / [2√π·(α-1)·K̄·sqrt(B)] · (1-Y²),  Y² = B/A
    Real Theta2 = (K_1_bar + M*mu_bar + kappa_bar) / (K_bar + M*mu_bar + kappa_bar);
    Real Hcoeff = (1.0 - Theta2) / (2.0*std::sqrt(M_PI)*(alpha-1.0)*K_bar*std::sqrt(B));
    epr = std::cos(-4.0*k_min_bar);
    epi = std::sin(-4.0*k_min_bar);
    Real m1pr = 1.0, m1pi = 1.0;
    Real Hr = m1pr*epr - m1pi*epi;
    Real Hi = m1pr*epi + m1pi*epr;
    Hr *= Hcoeff; Hi *= Hcoeff;

    // f2 = H·(Ẽ - e^{2iD} + i·(D + K̄ + Mμ̄ - κ̄)·G)
    epr = std::cos(2.0*D);
    epi = std::sin(2.0*D);
    Real part1r = Efr - epr;
    Real part1i = Efi - epi;

    Real coeffI = D + K_bar + M*mu_bar - kappa_bar;
    Real part2r = -coeffI*G_i;
    Real part2i =  coeffI*G_r;

    Real totalr = part1r + part2r;
    Real totali = part1i + part2i;

    f2r = Hr*totalr - Hi*totali;
    f2i = Hr*totali + Hi*totalr;
}


// R&M Eq. 16 — subcritical back-scattering correction.
//
//   A′₁ = K̄₁ + Mμ̄ − iκ̄′,   A′ = K̄ + Mμ̄ − iκ̄′,   Y′ = sqrt(A′₁/A′)
//   H′  = (1+i)(1 − Y′²) / (2√π·(α−1)·K̄·sqrt(A′₁))
//
//   I₂′ = (e^{−2iA′₁}/A′₁)·H′·{ A′·( e^{2iA′₁}·[1 − erf(sqrt(4κ̄′))] − 1 )
//                               + sqrt(2κ̄′)·(K̄ + Mμ̄ − μ̄·x₁/S0)
//                                 · F*(sqrt(−2iA′₁)) / sqrt(−iA′₁) }
//
// with F*(sqrt(i·w)) = (1+i)·E*(w), so F*(sqrt(−2iA′₁)) = (1+i)·E*(−2A′₁).
// All complex square roots are principal branch (complex_sqrt), matching the
// supercritical code.
//
// SIGNS ADJUDICATED NUMERICALLY, not read off (see CHANGELOG for the evidence).
// The exponent signs in the prefactor and in the erf bracket are inverted
// relative to the transcription this was implemented from: with e^{+2iA′₁}
// outside, Im(A′₁) = −κ̄′ makes the prefactor grow as e^{+2κ̄′} and |I| diverges
// (measured 4×10³⁷ at ξ = 20 instead of decaying). With the signs as written
// here every term decays, |I| is monotone decreasing in ξ, and the decrease is
// steeper at the higher frequency — the paper's Fig. 11 invariants.
template<typename Real>
void Radiation_integral2_subcrit(
    Real kappa_p, Real K_bar, Real mu_bar, Real S0,
    Real K_1_bar, Real alpha, Real x, Real M,
    Real &f2r, Real &f2i)
{
    Real A1p_r = K_1_bar + M*mu_bar,  A1p_i = -kappa_p;
    Real Ap_r  = K_bar   + M*mu_bar,  Ap_i  = -kappa_p;

    // Y′² = A′₁/A′
    Real Y2r, Y2i;
    cdiv<Real>(A1p_r, A1p_i, Ap_r, Ap_i, Y2r, Y2i);

    // H′ = (1+i)(1 − Y′²) / (2√π (α−1) K̄ sqrt(A′₁))
    Real om_r = 1.0 - Y2r, om_i = -Y2i;
    Real num_r = om_r - om_i;          // (1+i)·(1−Y′²)
    Real num_i = om_r + om_i;
    Real sA1_r, sA1_i; complex_sqrt<Real>(A1p_r, A1p_i, sA1_r, sA1_i);
    Real dscale = 2.0*std::sqrt(M_PI)*(alpha - 1.0)*K_bar;
    Real Hr, Hi;
    cdiv<Real>(num_r, num_i, dscale*sA1_r, dscale*sA1_i, Hr, Hi);

    // e^{−2iA′₁} = e^{−2κ̄′}·e^{−2i(K̄₁+Mμ̄)}, since Im(A′₁) = −κ̄′.
    Real decay = std::exp(-2.0*kappa_p);
    Real ex_r = decay*std::cos(-2.0*A1p_r);
    Real ex_i = decay*std::sin(-2.0*A1p_r);

    // term1 = (e^{−2iA′₁}/A′₁)·A′·( e^{2iA′₁}·[1 − erf(sqrt(4κ̄′))] − 1 )
    // collapsed analytically to  (A′/A′₁)·( [1 − erf(sqrt(4κ̄′))] − e^{−2iA′₁} ),
    // which is algebraically identical but never forms e^{+2κ̄′}. Evaluating the
    // two exponentials separately would overflow×underflow to NaN at large κ̄′.
    Real erf_r, erf_i;
    errFunc<Real>(std::sqrt(4.0*kappa_p), Real(0.0), erf_r, erf_i);
    Real b1_r = (1.0 - erf_r) - ex_r;
    Real b1_i = (    - erf_i) - ex_i;
    Real rat_r, rat_i;
    cdiv<Real>(Ap_r, Ap_i, A1p_r, A1p_i, rat_r, rat_i);      // A′/A′₁
    Real t1_r = rat_r*b1_r - rat_i*b1_i;
    Real t1_i = rat_r*b1_i + rat_i*b1_r;

    // term2 = (e^{−2iA′₁}/A′₁)·sqrt(2κ̄′)·(K̄ + Mμ̄ − μ̄x₁/S0)
    //         · F*(sqrt(−2iA′₁)) / sqrt(−iA′₁)
    //
    // F*(sqrt(−2iA′₁)) = (1+i)·E*(−2A′₁) = erf(ζ),  ζ = (1+i)·sqrt(−A′₁),
    // and ζ² = 2i·(−A′₁) = −2iA′₁ exactly, so the leading e^{−2iA′₁} is e^{ζ²}.
    // erf(ζ) grows as e^{+2κ̄′} and e^{ζ²} decays as e^{−2κ̄′}: formed separately
    // the product overflows to NaN for κ̄′ ≳ 350 (reachable from the validation
    // entry point at large chord × high frequency). Collapse it exactly instead:
    //     e^{ζ²}·erf(ζ) = e^{ζ²}(1 − erfc(ζ)) = e^{ζ²} − erfcx(ζ),
    // where both terms are bounded (erfcx(ζ) → −i/(Im ζ·√π) along the κ̄′ → ∞
    // asymptote ζ → i·sqrt(2κ̄′)).
    Real zr, zi;
    {   Real sr, si; complex_sqrt<Real>(-A1p_r, -A1p_i, sr, si);  // sqrt(−A′₁)
        zr = sr - si;                                             // (1+i)·sqrt(−A′₁)
        zi = sr + si;
    }
    Real ecx_r, ecx_i; erfcxFunc<Real>(zr, zi, ecx_r, ecx_i);
    Real W_r = ex_r - ecx_r;      // e^{−2iA′₁}·F*(sqrt(−2iA′₁))
    Real W_i = ex_i - ecx_i;

    Real sm_r, sm_i; complex_sqrt<Real>(A1p_i, -A1p_r, sm_r, sm_i);  // sqrt(−i·A′₁)
    Real den_r = A1p_r*sm_r - A1p_i*sm_i;      // A′₁·sqrt(−iA′₁)
    Real den_i = A1p_r*sm_i + A1p_i*sm_r;
    Real q_r, q_i; cdiv<Real>(W_r, W_i, den_r, den_i, q_r, q_i);
    Real c2 = std::sqrt(2.0*kappa_p) * (K_bar + M*mu_bar - mu_bar*x/S0);
    Real t2_r = c2*q_r, t2_i = c2*q_i;

    // I₂′ = H′ · (term1 + term2)
    Real s_r = t1_r + t2_r, s_i = t1_i + t2_i;
    f2r = Hr*s_r - Hi*s_i;
    f2i = Hr*s_i + Hi*s_r;
}


// |I| on whichever branch ξ selects, with no near-cutoff bridging.
// ξ < 1 supercritical (κ̄ = μ̄·sqrt(1−ξ²)); ξ > 1 subcritical (κ̄′ = μ̄·sqrt(ξ²−1)).
template<typename Real>
Real Amiet_I_abs_raw(Real xi, Real mu_bar, Real K_bar, Real K_1_bar, Real C,
                     Real S0, Real x, Real M, Real alpha)
{
    Real Ir, Ii;
    if (xi < 1.0) {
        Real kappa_bar = mu_bar * std::sqrt(1.0 - xi*xi);
        Real B = K_1_bar + M*mu_bar + kappa_bar;
        Real f1r, f1i;
        Radiation_integral1_cplx<Real>(B, Real(0.0), C, f1r, f1i);
        Real f2r, f2i;
        Radiation_integral2_general<Real>(B, K_bar, kappa_bar, mu_bar, S0,
                                          K_1_bar, alpha, x, M, f2r, f2i);
        Ir = f1r + f2r; Ii = f1i + f2i;
    } else {
        Real kappa_p = mu_bar * std::sqrt(xi*xi - 1.0);
        // Eq. 15 = Eq. 13 continued: B -> A′₁ = K̄₁ + Mμ̄ − iκ̄′.
        Real f1r, f1i;
        Radiation_integral1_cplx<Real>(K_1_bar + M*mu_bar, -kappa_p, C, f1r, f1i);
        Real f2r, f2i;
        Radiation_integral2_subcrit<Real>(kappa_p, K_bar, mu_bar, S0,
                                          K_1_bar, alpha, x, M, f2r, f2i);
        Ir = f1r + f2r; Ii = f1i + f2i;
    }
    return std::sqrt(Ir*Ir + Ii*Ii);
}


// |I| with the near-cutoff regularisation of R&M §4.1.
//
// The paper states only the procedure — "matching the values of the derivative
// ∂I/∂K₂ from both sides of the cuts and then re-calculating I" — and specifies
// no window or interpolant. THIS IS A CONCRETE REALISATION consistent with it,
// not a verbatim transcription:
//   * work in ξ; anchor at ξ_a = sqrt(1 − (κ̄_reg/μ̄)²) on the supercritical side
//     and ξ_b = sqrt(1 + (κ̄_reg/μ̄)²) on the subcritical side;
//   * at each anchor take |I| from its own branch and d|I|/dξ by a one-sided
//     finite difference stepping AWAY from the cut;
//   * cubic Hermite in ξ on |I| between them. On |I|, not on complex I: only
//     |I|² enters Eq. 18, and interpolating phase across the cut is meaningless.
//
// Low-frequency degeneracy: when μ̄ ≤ 2κ̄_reg the window would swallow ξ = 0, so
// κ̄_reg is capped at μ̄/2, keeping ξ_a ≥ sqrt(3)/2 and the anchors distinct.
//
// In production ξ = β|x₂|/S0 < 1, so only the ξ_a < ξ < 1 half of the bridge is
// ever reached; the subcritical anchor at ξ_b is what makes it well-posed.
template<typename Real>
Real Amiet_I_abs(Real xi, Real mu_bar, Real K_bar, Real K_1_bar, Real C,
                 Real S0, Real x, Real M, Real alpha)
{
    Real kr = Real(AMIET_KAPPA_REG);
    if (mu_bar <= 2.0*AMIET_KAPPA_REG) kr = 0.5*mu_bar;
    Real t = kr / mu_bar;                       // ∈ (0, 0.5]
    Real xi_a = std::sqrt(1.0 - t*t);
    Real xi_b = std::sqrt(1.0 + t*t);

    if (xi <= xi_a || xi >= xi_b)
        return Amiet_I_abs_raw<Real>(xi, mu_bar, K_bar, K_1_bar, C, S0, x, M, alpha);

    Real L = xi_b - xi_a;
    Real h = 0.01 * L;

    Real Ia   = Amiet_I_abs_raw<Real>(xi_a,     mu_bar, K_bar, K_1_bar, C, S0, x, M, alpha);
    Real Iam  = Amiet_I_abs_raw<Real>(xi_a - h, mu_bar, K_bar, K_1_bar, C, S0, x, M, alpha);
    Real Ib   = Amiet_I_abs_raw<Real>(xi_b,     mu_bar, K_bar, K_1_bar, C, S0, x, M, alpha);
    Real Ibp  = Amiet_I_abs_raw<Real>(xi_b + h, mu_bar, K_bar, K_1_bar, C, S0, x, M, alpha);

    Real da = (Ia - Iam) / h;     // one-sided, stepping away from the cut
    Real db = (Ibp - Ib) / h;

    Real s  = (xi - xi_a) / L;
    Real s2 = s*s, s3 = s2*s;
    Real h00 =  2.0*s3 - 3.0*s2 + 1.0;
    Real h10 =      s3 - 2.0*s2 + s;
    Real h01 = -2.0*s3 + 3.0*s2;
    Real h11 =      s3 -     s2;
    return h00*Ia + h10*L*da + h01*Ib + h11*L*db;
}


template<typename Real>
void Radiation_integral_total_vec(
    const std::vector<Real>& C,
    const std::vector<Real>& K_bar,
    const std::vector<Real>& mu_bar,
    Real S0,
    const std::vector<Real>& K_1_bar,
    const std::vector<Real>& K_2_bar,
    Real beta,
    Real alpha,
    Real M,
    Real x,
    std::vector<Real>& I_abs2
)
{
    for (std::size_t i = 0; i < C.size(); ++i)
    {
        // ξ = |K̄₂|/(βμ̄) = β|x₂|/S0. Only ξ² matters downstream, so take |K̄₂|.
        Real xi = std::abs(K_2_bar[i]) / (beta * mu_bar[i]);
        Real Iabs = Amiet_I_abs<Real>(xi, mu_bar[i], K_bar[i], K_1_bar[i], C[i],
                                      S0, x, M, alpha);
        I_abs2[i] = Iabs*Iabs;
    }
}


template<typename Real>
void TE_noise_outer_vec(
    Real M, Real U, Real x, Real y, Real z,
    Real b, Real c, Real span,

    Real c0,

    const std::vector<Real>& omega,
    const Real Ue_bot, const Real Ue_top,
    std::vector<Real>& WPS_lower, std::vector<Real>& WPS_upper,

    std::vector<Real>& farfieldSpectra
)
{
    const std::size_t N = omega.size();
    Real Ue[2] = {Ue_top, Ue_bot};

    Real beta = std::sqrt(1.0 - M*M);
    // General convected distance S0² = x₁² + β²(x₂² + x₃²). The β²y² term is
    // appended (rather than folded into a β²(y²+z²) group) so that at y = 0 it
    // is exactly +0.0 and the leading expression is textually the mid-span one
    // — making the mid-span reduction bit-identical here rather than merely
    // accurate to round-off.
    Real S0   = std::sqrt(x*x + beta*beta*z*z + beta*beta*y*y);

    for (int surf = 0; surf < 2; ++surf) {

    Real U_c  = 0.7 * Ue[surf];
    Real alpha = U / U_c;

    std::vector<Real> C(N);
    std::vector<Real> K_bar(N);
    std::vector<Real> mu_bar(N);
    std::vector<Real> K_1_bar(N);
    std::vector<Real> K_2_bar(N);
    std::vector<Real> K_2(N);       // dimensional spanwise wavenumber, for l_y

    for (std::size_t i = 0; i < N; ++i)
    {
        Real K = omega[i] / U;
        Real k = omega[i] / c0;     // acoustic wavenumber

        K_bar[i]   = K * b;
        mu_bar[i]  = K_bar[i] * M / (beta*beta);
        K_1_bar[i] = alpha * K_bar[i];
        C[i]       = K_1_bar[i] - mu_bar[i] * (x/S0 - M);   // Eq. 12
        // Eq. 18's large-aspect-ratio delta selects the single gust K₂ = k·x₂/S0.
        K_2[i]     = k * y / S0;
        K_2_bar[i] = K_2[i] * b;
    }

    std::vector<Real> I_abs2(N);
    Radiation_integral_total_vec(C, K_bar, mu_bar, S0, K_1_bar, K_2_bar, beta,
                                 alpha, M, x, I_abs2);

    // Roger & Moreau (2005) Eq. 19 — spanwise-wavenumber-corrected Corcos length.
    Real b_c = 1.47;
    std::vector<Real> l_y(N);
    for (std::size_t i = 0; i < N; ++i)
        l_y[i] = corcos_l_y<Real>(K_2[i], omega[i], b_c, U_c);

    Real b_half = b;  // b is already the semi-chord (c/2) passed by caller
    for (std::size_t i = 0; i < N; ++i) {
        // Eq. 18's prefactor keeps the plate-normal x₃ in the numerator; only
        // S0 generalises.
        Real _t1   = (omega[i]*b_half*z) / (2.0*M_PI*c0*S0*S0);
        Real term1 = _t1 * _t1;

        if (surf == 0) {
            farfieldSpectra[i]  = term1 * 2.0*span * I_abs2[i] * WPS_upper[i] * l_y[i];
        } else {
            farfieldSpectra[i] += term1 * 2.0*span * I_abs2[i] * WPS_lower[i] * l_y[i];
        }
    }
    }

}
