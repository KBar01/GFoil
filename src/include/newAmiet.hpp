#pragma once

// Mid-span (x₂ = 0) trailing-edge noise model — Roger & Moreau (2005).
// Restriction: observer in the mid-span plane only.
//   x₂ = 0  ⟹  K̄₂ = 0  ⟹  κ̄ = μ̄
//   kbar² = μ̄² - (K̄₂/β)² = μ̄² > 0 always  ⟹  supercritical path only.
// The subcritical branch and all associated helpers are not present.

#include <cmath>
#include <complex>
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


/////////////////////////////////////////// Amiet model — Roger & Moreau (2005) ///////////////////

// Computes E*(x) used in Amiet transfer function, Roger & Moreau (2005) Eq.6.
// E(x) = conj(E*(x)) obtained by negating imaginary part.
// E*(xr + i·xi) = erf[(1+i)·sqrt((xr+i·xi)/2)] / (1+i)
// For real x (xi=0) this gives the mid-span E* of R&M Eq. 3-4.
template<typename Real>
inline void Fresnel_int_conj(Real xr, Real xi,
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
    Real a_r, a_i; Fresnel_int_conj<Real>(2.0*(B-C), 0.0, a_r, a_i);
    Real b_r, b_i; Fresnel_int_conj<Real>(2.0*B,     0.0, b_r, b_i);

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

    Real t1r = tmp_r*a_div_r - tmp_i*a_div_i;
    Real t1i = tmp_r*a_div_i + tmp_i*a_div_r;

    Real t2r = -(onepI_r*b_r - onepI_i*b_i);
    Real t2i = -(onepI_r*b_i + onepI_i*b_r);

    Real br_r = t1r + t2r + 1.0;
    Real br_i = t1i + t2i;

    f1r = pref_r*br_r - pref_i*br_i;
    f1i = pref_r*br_i + pref_i*br_r;
}


// Back-scattering correction, Roger & Moreau (2005) Eq.14.
// G is the sum of sub-integrals G_a..G_e; H is the correction prefactor; ε from Eq.9.
template<typename Real>
void Radiation_integral2(
    Real B, Real K_bar, Real k_min_bar, Real mu_bar, Real S0,
    Real K_1_bar, Real alpha, Real x, Real M,
    Real &f2r, Real &f2i)
{
    Real error = std::pow(1.0 + 1.0/(4.0*mu_bar), -0.5);   // Eq.9
    // D = mu_bar*(1 - x/S0); k_min_bar = mu_bar at mid-span (K_2_bar=0)
    Real D = mu_bar * (1.0 - x / S0);

    // Ẽ = exp(4i·k_min_bar)·(1 - (1+i)·E*(4·k_min_bar))
    Real Fr, Fi;
    Fresnel_int_conj<Real>(4.0*k_min_bar, 0.0, Fr, Fi);
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
    Real denC = 2.0 * ((std::abs(denC_val) < Real(1e-10)) ? Real(1e-10) : denC_val);
    Real m1r = 1.0, m1i = -1.0;  // (1-i)
    Real coeffr = (1.0+error)*m1r / denC;
    Real coeffi = (1.0+error)*m1i / denC;
    epr = std::cos(4.0*k_min_bar);
    epi = std::sin(4.0*k_min_bar);
    Fresnel_int_conj<Real>(4.0*k_min_bar, 0.0, Fr, Fi);
    Real tmp_r = epr*Fr - epi*Fi;
    Real tmp_i = epr*Fi + epi*Fr;
    Real G_cr = coeffr*tmp_r - coeffi*tmp_i;
    Real G_ci = coeffr*tmp_i + coeffi*tmp_r;

    // --- G_d: [(1-ε)(1+i)] / [2(D+2k)] · e^{-4ik}·E(4k),  subtracted in sum
    Real denD_val = D + 2.0*k_min_bar;
    Real denD = 2.0 * ((std::abs(denD_val) < Real(1e-10)) ? Real(1e-10) : denD_val);
    Real p1r = 1.0, p1i = 1.0;  // (1+i)
    Real coeffDr = (1.0-error)*p1r / denD;
    Real coeffDi = (1.0-error)*p1i / denD;
    epr = std::cos(-4.0*k_min_bar);
    epi = std::sin(-4.0*k_min_bar);
    // E(4k) = conj(E*(4k)): compute E* then negate imaginary part
    Fresnel_int_conj<Real>(4.0*k_min_bar, 0.0, Fr, Fi);
    Fi = -Fi;
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
        Real sqrtfactor_r, sqrtfactor_i;
        complex_sqrt<Real>(k_min_bar / D, Real(0.0), sqrtfactor_r, sqrtfactor_i);

        Fresnel_int_conj<Real>(2.0*D, 0.0, Fr, Fi);
        tmp_r = e2r*Fr - e2i*Fi;
        tmp_i = e2r*Fi + e2i*Fr;
        Real new_r = tmp_r*sqrtfactor_r - tmp_i*sqrtfactor_i;
        Real new_i = tmp_r*sqrtfactor_i + tmp_i*sqrtfactor_r;
        tmp_r = new_r; tmp_i = new_i;

        // bracket = (1+i)·(1-ε)/(D+2k) - (1-i)·(1+ε)/(D-2k)
        Real da = D + 2.0*k_min_bar;
        Real db = D - 2.0*k_min_bar;
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



// Top-level Amiet TE noise. Computes far-field PSD at observer (x,y,z)
// for each frequency in omega[Nsound]. Roger & Moreau (2005) Eqs.1-2, 18.
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

    for (int surf = 0; surf < 2; ++surf) {

    Real beta = std::sqrt(1.0 - M*M);
    Real S0   = std::sqrt(x*x + beta*beta*(y*y + z*z));

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
    // S_pp = (ωb·z / 2πc₀S₀²)² · 2·span · |I|² · Φ_pp · l_y,  b = c/2 (half-chord).
    Real b_half = c / 2.0;
    for (int i = 0; i < Nsound; ++i) {
        Real term1 = std::pow((omega[i]*b_half*z) / (2.0*M_PI*c0*S0*S0), 2.0);

        if (surf == 0) {
            farfieldSpectra[i]  = term1 * 2.0*span * I_abs2[i] * WPS_upper[i] * l_y[i];
        } else {
            farfieldSpectra[i] += term1 * 2.0*span * I_abs2[i] * WPS_lower[i] * l_y[i];
        }
    }
    }

}
