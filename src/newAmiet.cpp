#include <cmath>
#include <complex>
#include "Faddeeva.hh"
#include "real_type.h"


using std::complex;
using std::exp;

// multiply (ar + i ai)*(br + i bi)
inline void cmul(Real ar, Real ai, Real br, Real bi, Real &cr, Real &ci) {
    cr = ar*br - ai*bi;
    ci = ar*bi + ai*br;
}

// divide (ar + i ai)/(br + i bi)
inline void cdiv(Real ar, Real ai, Real br, Real bi, Real &cr, Real &ci) {
    Real den = br*br + bi*bi;
    cr = (ar*br + ai*bi)/den;
    ci = (ai*br - ar*bi)/den;
}

// exp(ar+i ai)
inline void cexp(Real ar, Real ai, Real &er, Real &ei) {
    Real e = std::exp(ar);
    er = e*std::cos(ai);
    ei = e*std::sin(ai);
}

// sqrt(ar+i ai)
inline void csqrt(Real ar, Real ai, Real &sr, Real &si) {
    Real r = std::sqrt(std::hypot(ar,ai));
    Real t = std::atan2(ai,ar)/2.0;
    sr = std::sqrt(r)*std::cos(t);
    si = std::sqrt(r)*std::sin(t);
}


inline void complex_sqrt(Real ar, Real ai,
                         Real &br, Real &bi)
{
    // Computes sqrt(ar + i ai) = br + i bi
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

    Real r = std::hypot(ar, ai);   // sqrt(ar^2 + ai^2)
    Real t = std::sqrt(0.5*(r + std::fabs(ar)));

    if (ar >= 0.0) {
        br = t;
        bi = ai/(2.0*t);
    } else {
        br = std::fabs(ai)/(2.0*t);
        bi = (ai >= 0.0) ? t : -t;
    }
}


void errFunc(Real in_r,Real in_i, Real &out_r, Real &out_i){

    double x_val = in_r.getValue() ;
    //double x_val = in_r ;
    double y_val = in_i.getValue() ;
    //double y_val = in_i ;
    complex<double> z(x_val,y_val) ;
    complex<double> w = Faddeeva::erf(z);
    
    std::complex<double> dw_dz = (2.0 / sqrt(M_PI)) * exp(-z * z); // Derivative

    // Compute real and imag parts of Jacobian
    double du_dx = dw_dz.real();        // ∂Re(w)/∂x
    double du_dy = -dw_dz.imag();       // ∂Re(w)/∂y
    double dv_dx = dw_dz.imag();        // ∂Im(w)/∂x
    double dv_dy = dw_dz.real();        // ∂Im(w)/∂y
    
    // Push statement for u = Re(w)
    codi::StatementPushHelper<Real> ph;
    ph.startPushStatement();
    ph.pushArgument(in_r, du_dx);
    ph.pushArgument(in_i, du_dy);
    ph.endPushStatement(out_r, w.real());


    // Push statement for v = Im(w)

    codi::StatementPushHelper<Real> phIm;
    phIm.startPushStatement();
    phIm.pushArgument(in_r, dv_dx);
    phIm.pushArgument(in_i, dv_dy);
    phIm.endPushStatement(out_i, w.imag());
}


/////////////////////////////////////////// This is amiet model as given in R&M ///////////////////

void wavesnumbers(
    // wavesnumbers: compute wavenumber-related quantities from input scalars

    Real M, Real U, Real x, Real y, Real z,  // mach,freestream, observer x,y,x from TE
    Real b, Real Ky,
    Real c0,
    // Frequency array (radians/sec):
    const Real omega[Nsound],

    // Outputs (arrays):
    Real C[Nsound], // constant
    Real K_bar[Nsound],
    Real mu_bar[Nsound],
    Real K_1_bar[Nsound],
    Real K_2_bar[Nsound],
    Real kappa_bar[Nsound],

    // Outputs (scalars):
    Real S0,
    Real alpha,
    Real U_c,
    Real beta
)
{
    // acoustic wavenumber array
    // kappa_bar will be filled inside the loop
    beta = std::sqrt(1.0 - M*M);

    // distance to observer
    S0 = std::sqrt(x*x + beta*beta * (y*y + z*z));

    // convection velocity
    U_c = 0.7 * U; 
    alpha = U / U_c;

    for (int i=0; i<Nsound; ++i)
    {
        // acoustic wavenumber
        Real kappa = omega[i] / c0;

        // aerodynamic wavenumber
        Real K = omega[i] / U;

        // non-dimensional forms
        K_bar[i]    = K * b;
        kappa_bar[i]= kappa * b;

        // mu_bar
        mu_bar[i]   = K_bar[i] * M / (beta*beta);

        // K_1_bar
        K_1_bar[i]  = alpha * K_bar[i];

        // C (temporal variable)
        C[i] = K_1_bar[i] - mu_bar[i] * (x / S0 - M);       // EQ 12

        // K_2_bar
        if (Ky == 0.0) {
            K_2_bar[i] = kappa_bar[i] * y / S0;
        } else {
            // If you really want the linspace(-Ky,Ky,100) branch,
            // you’d handle it separately. For now we mimic Ky==0 case.
            K_2_bar[i] = kappa_bar[i] * y / S0;
        }
    }
}


// x is real Real
void Fresnel_int(Real x,
                 Real &E_real,
                 Real &E_imag)
{
    // Step 1: complex argument z = (1 - i)*sqrt(x/2)
    Real s = std::sqrt(0.5 * x);
    Real z_real =  s;       // real part (1 - i)*s
    Real z_imag = -s;       // imag part

    // Step 2: call your complex erf wrapper
    Real erf_real, erf_imag;
    errFunc(z_real, z_imag, erf_real, erf_imag);

    // Step 3: divide by (1 - i)
    // (a+ib)/(1-i) = [(a+ib)(1+i)]/2
    Real denom = 2.0; // |1-i|^2
    E_real = ( erf_real*1.0 - erf_imag*1.0 ) / denom;
    E_imag = ( erf_real*1.0 + erf_imag*1.0 ) / denom;
}

// x is real (Real)
inline void Fresnel_int_conj(Real xr, Real xi,
                             Real &Er, Real &Ei)
{
    // input x = xr + i xi
    // output E_s = Er + i Ei

    // compute (1+1i)*sqrt(0.5*x)
    Real sr, si; 
    complex_sqrt(0.5*xr, 0.5*xi, sr, si); // reuse the sqrt we wrote earlier
    // multiply (1+1i)*(sr+isi)
    Real mr = (1.0)*sr - (1.0)*si; // real part
    Real mi = (1.0)*sr + (1.0)*si; // imag part

    // erfz of that complex number
    Real erf_r, erf_i;
    errFunc(mr, mi, erf_r, erf_i);  // you provide this

    // temp = -(1 - erfz(...))
    Real temp_r = -(1.0 - erf_r);
    Real temp_i = -(-erf_i); // = +erf_i

    // temp + 1
    temp_r += 1.0;

    // divide by (1+1i): (a+ib)/(1+i) = ((a+b) + i(b-a))/2
    Real a = temp_r;
    Real b = temp_i;
    Er = 0.5*(a + b);
    Ei = 0.5*(b - a);
}

inline void Phi_0_img_new(Real sqrt_r, Real sqrt_i,
                          Real &Phi_r, Real &Phi_i)
{
    // input: sqrt_r + i sqrt_i = sqrt(i*x)
    // output: Phi_0 = Phi_r + i Phi_i

    // x = sqrt_ix^2 * (-i)
    Real a = sqrt_r*sqrt_r - sqrt_i*sqrt_i;     // real part of sqrt_ix^2
    Real b = 2.0*sqrt_r*sqrt_i;                 // imag part of sqrt_ix^2
    // multiply by -i: (a+ib)*(-i) = b - i a
    Real xr = b;
    Real xi = -a;

    // Fresnel_int_conj(x) returning Er+ iEi
    Real Er, Ei;
    Fresnel_int_conj(xr, xi, Er, Ei);

    // (1+i)*(Er+iEi) = (Er - Ei) + i(Er+Ei)
    Phi_r = Er - Ei;
    Phi_i = Er + Ei;
}

inline void Radiation_integral1(Real B, Real C,
                                Real &f1r, Real &f1i)
{
    // Fresnel integrals
    Real a_r,a_i; Fresnel_int_conj(2.0*(B-C),0.0,a_r,a_i);
    Real b_r,b_i; Fresnel_int_conj(2.0*B,0.0,b_r,b_i);

    // prefactor = -exp(2 i C)/(i C)
    Real cos2C=std::cos(2.0*C), sin2C=std::sin(2.0*C);
    // exp(2iC)=cos2C+i sin2C
    // (iC)=i*C so 1/(iC)=-i/C
    // so prefactor = -exp(2iC)/(iC)= -exp(2iC)*(-i/C)= i*exp(2iC)/C
    Real pref_r = -sin2C/C;   // real part of i*exp(2iC)/C
    Real pref_i =  cos2C/C;   // imag part

    // (1+i)
    Real onepI_r=1.0, onepI_i=1.0;

    // exp(-2iC)
    Real e_2C_r=std::cos(-2.0*C), e_2C_i=std::sin(-2.0*C); // =cos2C - i sin2C

    // sqrt(2B)
    Real s = std::sqrt(2.0*B);

    // term1 = (1+i)*exp(-2iC)*s*a/sqrt(2(B-C))
    //Real sc=std::sqrt(2.0*(B-C));
    Real sc_r, sc_i;
    complex_sqrt(2.0*(B-C), 0.0, sc_r, sc_i);

    // first multiply (1+i)*exp(-2iC)
    Real tmp_r=onepI_r*e_2C_r - onepI_i*e_2C_i;
    Real tmp_i=onepI_r*e_2C_i + onepI_i*e_2C_r;
    // multiply by s
    tmp_r*=s; tmp_i*=s;
    // divide a/sc
    //Real a_div_r=a_r/sc, a_div_i=a_i/sc;
    Real a_div_r, a_div_i;
    cdiv(a_r, a_i, sc_r, sc_i, a_div_r, a_div_i);
    
    
    // multiply tmp * a_div
    Real t1r=tmp_r*a_div_r - tmp_i*a_div_i;
    Real t1i=tmp_r*a_div_i + tmp_i*a_div_r;

    // term2 = -(1+i)*b
    Real t2r=-(onepI_r*b_r - onepI_i*b_i);
    Real t2i=-(onepI_r*b_i + onepI_i*b_r);

    // bracket = term1+term2+1
    Real br_r=t1r+t2r+1.0;
    Real br_i=t1i+t2i;

    // f1 = prefactor*bracket
    f1r=pref_r*br_r - pref_i*br_i;
    f1i=pref_r*br_i + pref_i*br_r;
}


void Radiation_integral2(
    Real B, Real K_bar, Real k_min_bar, Real mu_bar, Real S0,
    Real K_1_bar, Real alpha, Real x, Real M,
    /*outputs:*/ Real &f2r, Real &f2i)
{
    Real error = std::pow(1.0 + 1.0/(4.0*mu_bar), -0.5);   // Eq.9
    Real D = k_min_bar - mu_bar * x / S0;

    // --- E = exp(4i*k_min_bar)*(1 - (1+i)*Fresnel_int_conj(4*k_min_bar))
    Real Fr,Fi;
    Fresnel_int_conj(4.0*k_min_bar,0.0,Fr,Fi);
    // (1+i)*Fresnel
    Real t1r=1.0*Fr -1.0*Fi;
    Real t1i=1.0*Fi +1.0*Fr;
    // 1 - that
    Real oneMinus_r=1.0 - t1r;
    Real oneMinus_i=0.0 - t1i;
    // exp(4i*k_min_bar)
    Real e4r=std::cos(4.0*k_min_bar);
    Real e4i=std::sin(4.0*k_min_bar);
    Real Er=e4r*oneMinus_r - e4i*oneMinus_i;
    Real Ei=e4r*oneMinus_i + e4i*oneMinus_r;

    // E_final = Re(E) + i*error*Im(E)
    Real Efr=Er;
    Real Efi=error*Ei;

    // --- G_a
    Real phase=2.0*k_min_bar+D;
    Real epr=std::cos(phase);
    Real epi=std::sin(phase);
    Real G_ar=(1.0+error)*epr*std::sin(D-2.0*k_min_bar)/(D-2.0*k_min_bar);
    Real G_ai=(1.0+error)*epi*std::sin(D-2.0*k_min_bar)/(D-2.0*k_min_bar);

    // --- G_b
    phase=-2.0*k_min_bar+D;
    epr=std::cos(phase);
    epi=std::sin(phase);
    Real G_br=(1.0-error)*epr*std::sin(D+2.0*k_min_bar)/(D+2.0*k_min_bar);
    Real G_bi=(1.0-error)*epi*std::sin(D+2.0*k_min_bar)/(D+2.0*k_min_bar);

    // --- G_c
    Real denC=2.0*(D-2.0*k_min_bar);
    // (1- i)=1 - i
    Real m1r=1.0, m1i=-1.0;
    Real coeffr=(1.0+error)*m1r/denC;
    Real coeffi=(1.0+error)*m1i/denC;
    epr=std::cos(4.0*k_min_bar);
    epi=std::sin(4.0*k_min_bar);
    Fresnel_int_conj(4.0*k_min_bar,0.0,Fr,Fi);
    // multiply exp * Fresnel
    Real tmp_r=epr*Fr - epi*Fi;
    Real tmp_i=epr*Fi + epi*Fr;
    // then *coeff
    Real G_cr=coeffr*tmp_r - coeffi*tmp_i;
    Real G_ci=coeffr*tmp_i + coeffi*tmp_r;

    // --- G_d
    Real denD=2.0*(D+2.0*k_min_bar);
    Real p1r=1.0, p1i=1.0;
    Real coeffDr=(1.0-error)*p1r/denD;
    Real coeffDi=(1.0-error)*p1i/denD;
    epr=std::cos(-4.0*k_min_bar);
    epi=std::sin(-4.0*k_min_bar);
    Fresnel_int(4.0*k_min_bar,Fr,Fi);
    tmp_r=epr*Fr - epi*Fi;
    tmp_i=epr*Fi + epi*Fr;
    Real G_dr=coeffDr*tmp_r - coeffDi*tmp_i;
    Real G_di=coeffDr*tmp_i + coeffDi*tmp_r;

    // --- G_e
    Real e2r=std::cos(2.0*D);
    Real e2i=std::sin(2.0*D);
    //Real sqrtfactor=std::sqrt(2.0*k_min_bar)/std::sqrt(2.0*D);
    // numerator: 2*k_min_bar (assume real)
    Real num_r = 2.0 * k_min_bar;
    Real num_i = 0.0;
    // denominator: 2*D (assume real here but allow 0 imag)
    Real den_r = 2.0 * D;
    Real den_i = 0.0;
    // divide num/den as complex:
    Real den_mag2 = den_r*den_r + den_i*den_i;
    Real quot_r = (num_r*den_r + num_i*den_i) / den_mag2;
    Real quot_i = (num_i*den_r - num_r*den_i) / den_mag2;

    // now sqrt that quotient:
    Real sqrtfactor_r, sqrtfactor_i;
    complex_sqrt(quot_r, quot_i, sqrtfactor_r, sqrtfactor_i);



    Fresnel_int_conj(2.0*D,0.0,Fr,Fi);
    // multiply e2 *Fresnel*sqrtfactor
    tmp_r=e2r*Fr - e2i*Fi;
    tmp_i=e2r*Fi + e2i*Fr;

    Real new_r = tmp_r*sqrtfactor_r - tmp_i*sqrtfactor_i;
    Real new_i = tmp_r*sqrtfactor_i + tmp_i*sqrtfactor_r;
    tmp_r=new_r; tmp_i=new_i;

    // inside big bracket:
    Real a1r=(1.0+1.0)*( (1.0-error)/(D+2.0*k_min_bar)); // (1+i) part is separate
    // Actually bracket = (1+i)*(...) - (1-i)*(...)
    Real term1r=(1.0)*((1.0-error)/(D+2.0*k_min_bar));
    Real term1i=(1.0)*((1.0-error)/(D+2.0*k_min_bar));
    Real term2r=(1.0)*((1.0+error)/(D-2.0*k_min_bar));
    Real term2i=-1.0*((1.0+error)/(D-2.0*k_min_bar));
    // (1+i)*A - (1-i)*B
    Real Br_r=term1r - term2r;
    Real Br_i=term1i - term2i;
    // multiply tmp * bracket
    Real G_er=tmp_r*Br_r - tmp_i*Br_i;
    Real G_ei=tmp_r*Br_i + tmp_i*Br_r;

    // sum
    Real G_r=G_ar+G_br+G_cr -G_dr+G_er;
    Real G_i=G_ai+G_bi+G_ci -G_di+G_ei;

    // --- H
    Real Theta=std::sqrt( (K_1_bar+(1.0+M)*mu_bar)/(K_bar+(1.0+M)*mu_bar));
    Real Hcoeff=(1.0-Theta*Theta)/(2.0*std::sqrt(M_PI)*(alpha-1.0)*K_bar*std::sqrt(B));
    epr=std::cos(-4.0*k_min_bar);
    epi=std::sin(-4.0*k_min_bar);
    Real m1pr=1.0, m1pi=1.0; // (1+i)
    Real Hr=m1pr*epr - m1pi*epi;
    Real Hi=m1pr*epi + m1pi*epr;
    Hr*=Hcoeff; Hi*=Hcoeff;

    // final f2 = H*(E_final - exp(2iD) + i*(D+K_bar+M*mu_bar-k_min_bar)*G)
    // part1 = E_final - exp(2iD)
    epr=std::cos(2.0*D);
    epi=std::sin(2.0*D);
    Real part1r=Efr - epr;
    Real part1i=Efi - epi;

    // part2 = i*(D+...)*G
    Real coeffI=D+K_bar+M*mu_bar-k_min_bar;
    // i*coeff *G = -coeff*G_i + i*coeff*G_r
    Real part2r=-coeffI*G_i;
    Real part2i= coeffI*G_r;

    Real totalr=part1r+part2r;
    Real totali=part1i+part2i;

    // multiply by H
    f2r=Hr*totalr - Hi*totali;
    f2i=Hr*totali + Hi*totalr;
}

// your main function:
void Radiation_integral1_subcrit(
    Real C,
    Real A1prime_r, Real A1prime_i,
    Real mu_bar,
    Real x,
    Real S0,
    Real k_min_bar_prime,
    Real &outReal,
    Real &outImag)
{
    // temp1 = 2*(mu_bar*(x/S0) - i*k_min_bar_prime)
    Real temp1_real = 2.0 * mu_bar * (x / S0);
    Real temp1_imag = -2.0 * k_min_bar_prime;

    // temp2 = 2*A1_prime (complex)
    Real temp2_r = 2.0 * A1prime_r;
    Real temp2_i = 2.0 * A1prime_i;

    // exp(±2iC)
    Real epos_real = std::cos(2.0 * C);
    Real epos_imag = std::sin(2.0 * C);
    Real eneg_real = std::cos(-2.0 * C);
    Real eneg_imag = std::sin(-2.0 * C);

    // 1/(i*C) = -i/C
    Real inv_iC_real = 0.0;
    Real inv_iC_imag = -1.0 / C;

    // (1+ i)
    Real onep_i_real = 1.0;
    Real onep_i_imag = 1.0;

    // Fresnel_int_conj(temp1)
    Real Fres_real, Fres_imag;
    Fresnel_int_conj(temp1_real, temp1_imag, Fres_real, Fres_imag);

    // sqrt(temp1) for denominator (complex)
    Real sqrt1_real, sqrt1_imag;
    complex_sqrt(temp1_real, temp1_imag, sqrt1_real, sqrt1_imag);

    // (1+i)*Fres
    Real t1_real = onep_i_real * Fres_real - onep_i_imag * Fres_imag;
    Real t1_imag = onep_i_real * Fres_imag + onep_i_imag * Fres_real;

    // divide by sqrt(temp1)
    Real denom_r = sqrt1_real;
    Real denom_i = sqrt1_imag;
    Real denom_mod2 = denom_r * denom_r + denom_i * denom_i;
    Real t1d_real = (t1_real * denom_r + t1_imag * denom_i) / denom_mod2;
    Real t1d_imag = (t1_imag * denom_r - t1_real * denom_i) / denom_mod2;

    // sqrt(temp2) (complex now)
    Real sqrt_temp2_r, sqrt_temp2_i;
    complex_sqrt(temp2_r, temp2_i, sqrt_temp2_r, sqrt_temp2_i);

    // multiply by sqrt(temp2)
    Real t1ds_r = t1d_real * sqrt_temp2_r - t1d_imag * sqrt_temp2_i;
    Real t1ds_i = t1d_real * sqrt_temp2_i + t1d_imag * sqrt_temp2_r;

    // multiply by exp(-2iC)
    Real partA_real = eneg_real * t1ds_r - eneg_imag * t1ds_i;
    Real partA_imag = eneg_real * t1ds_i + eneg_imag * t1ds_r;

    // phi_0_img_new( sqrt(2*A1_prime*i) )
    Real arg_r = -A1prime_i * 2.0; // real part of 2*A1_prime*i
    Real arg_i = A1prime_r * 2.0;  // imag part
    Real sqrt_arg_r, sqrt_arg_i;
    complex_sqrt(arg_r, arg_i, sqrt_arg_r, sqrt_arg_i);
    Real phi_real, phi_imag;
    Phi_0_img_new(sqrt_arg_r, sqrt_arg_i, phi_real, phi_imag);

    // inner = partA - phi + 1
    Real inner_real = partA_real - phi_real + 1.0;
    Real inner_imag = partA_imag - phi_imag;

    // multiply by (-exp(2iC)/(iC))
    Real e2_real = std::cos(2.0 * C);
    Real e2_imag = std::sin(2.0 * C);
    e2_real = -e2_real;
    e2_imag = -e2_imag;

    Real m_real = e2_real * inv_iC_real - e2_imag * inv_iC_imag;
    Real m_imag = e2_real * inv_iC_imag + e2_imag * inv_iC_real;

    outReal = m_real * inner_real - m_imag * inner_imag;
    outImag = m_real * inner_imag + m_imag * inner_real;
}

inline void Radiation_integral2_subcrit(
    Real A_prime_r, Real A_prime_i,
    Real A1_prime_r, Real A1_prime_i,
    Real mu_bar,
    Real M, Real x, Real S0,
    Real k_min_bar_prime,
    Real alpha, Real K_bar,
    Real Theta_prime_r, Real Theta_prime_i,
    Real &out_r, Real &out_i)
{
    // D' = mu_bar*x/S0 - i*k_min
    Real Dr = mu_bar*x/S0;
    Real Di = -k_min_bar_prime;

    // Theta^2 = (Θr+iΘi)^2
    Real Th2r = Theta_prime_r*Theta_prime_r - Theta_prime_i*Theta_prime_i;
    Real Th2i = 2.0*Theta_prime_r*Theta_prime_i;

    // (1 - Theta^2)
    Real OneMinus_r = 1.0 - Th2r;
    Real OneMinus_i = -Th2i;

    // (1+i)*(1-Theta^2)
    Real num_r = OneMinus_r - OneMinus_i;
    Real num_i = OneMinus_r + OneMinus_i;

    // sqrt(A1_prime) (complex)
    Real sqrtA1_r, sqrtA1_i;
    complex_sqrt(A1_prime_r,A1_prime_i,sqrtA1_r,sqrtA1_i);
    Real mod_sqrtA1 = std::sqrt(sqrtA1_r*sqrtA1_r + sqrtA1_i*sqrtA1_i);

    Real denom = 2.0*std::sqrt(M_PI)*(alpha-1.0)*K_bar*mod_sqrtA1;

    // H' = num/denom
    Real Hr = num_r/denom;
    Real Hi = num_i/denom;

    // exp(-2i*D') = exp(2Di - i2Dr)
    Real exp1_r = std::cos(-2.0*Dr)*std::exp(2.0*Di);
    Real exp1_i = std::sin(-2.0*Dr)*std::exp(2.0*Di);

    // 1/D'
    Real denom_sq = Dr*Dr+Di*Di;
    Real invDr =  Dr/denom_sq;
    Real invDi = -Di/denom_sq;

    // pre = exp(-2iD')/D'
    Real pre_r = exp1_r*invDr - exp1_i*invDi;
    Real pre_i = exp1_r*invDi + exp1_i*invDr;

    // exp(2i*D') = exp(-2Di + i2Dr)
    Real exp2_r = std::cos(2.0*Dr)*std::exp(-2.0*Di);
    Real exp2_i = std::sin(2.0*Dr)*std::exp(-2.0*Di);

    // erf(sqrt(4*k_min)) (real arg)
    Real sqarg = std::sqrt(4.0*k_min_bar_prime);
    Real erf_r,erf_i;
    errFunc(sqarg,0.0,erf_r,erf_i);

    // (1 - erf)
    Real omr = 1.0 - erf_r;
    Real omi = -erf_i;

    // part A = A_prime*exp(2iD')*(1 - erf)
    Real tmpAr = exp2_r*omr - exp2_i*omi;
    Real tmpAi = exp2_r*omi + exp2_i*omr;
    // multiply by A_prime (complex)
    Real newAr = A_prime_r*tmpAr - A_prime_i*tmpAi;
    Real newAi = A_prime_r*tmpAi + A_prime_i*tmpAr;
    tmpAr = newAr; tmpAi=newAi;

    // part B = 2*√(2*kmin)*(K_bar+(M-x/S0)*mu_bar)*ES_conj(temp3)
    Real temp3r = -2.0*Dr;
    Real temp3i =  2.0*Di;
    Real ESr,ESi;
    Fresnel_int_conj(temp3r,temp3i,ESr,ESi);

    Real coeff = 2.0*std::sqrt(2.0*k_min_bar_prime)*
                   (K_bar + (M-x/S0)*mu_bar);

    Real tmpBr = coeff*ESr;
    Real tmpBi = coeff*ESi;

    // inside brackets: partA -1 + partB
    Real br_r = tmpAr - 1.0 + tmpBr;
    Real br_i = tmpAi + tmpBi;

    // multiply by H'
    Real Hbr_r = Hr*br_r - Hi*br_i;
    Real Hbr_i = Hr*br_i + Hi*br_r;

    // multiply by pre
    out_r = pre_r*Hbr_r - pre_i*Hbr_i;
    out_i = pre_r*Hbr_i + pre_i*Hbr_r;
}

inline void Wavenumbers_subcrit(
    const Real M,
    Real mu_bar, Real K_1_bar, Real K_2_bar, Real K_bar, Real beta,
    Real &k_min_bar_prime,
    Real &A1prime_r, Real &A1prime_i,
    Real &Aprime_r, Real &Aprime_i,
    Real &Tr, Real &Ti)
{
    // k_min_bar_prime real
    Real arg = -mu_bar*mu_bar + (K_2_bar/beta)*(K_2_bar/beta);
    
    Real zero = 0.0;
    k_min_bar_prime = (arg > 0.0) ? std::sqrt(arg) : zero;

    // A1' = K1 + M*mu - i*kmin
    A1prime_r = K_1_bar + M*mu_bar;
    A1prime_i = -k_min_bar_prime;

    // A' = K + M*mu - i*kmin
    Aprime_r = K_bar + M*mu_bar;
    Aprime_i = -k_min_bar_prime;

    // ratio = A1'/A'
    Real denom = Aprime_r*Aprime_r + Aprime_i*Aprime_i;
    Real rr = (A1prime_r*Aprime_r + A1prime_i*Aprime_i)/denom;
    Real ri = (A1prime_i*Aprime_r - A1prime_r*Aprime_i)/denom;

    // Theta' = sqrt(rr + i ri)
    complex_sqrt(rr, ri, Tr, Ti);
}

void Radiation_integral_total(
    const Real *C,           // array size Nsound
    const Real *K_bar,       // array size Nsound
    const Real *mu_bar,      // array size Nsound
    Real S0,
    const Real *K_1_bar,     // array size Nsound
    const Real *K_2_bar,     // array size Nsound
    Real alpha,
    Real M,
    Real x,
    Real beta,
    Real Ky,
    Real *I_abs2             // output |I|^2, size Nsound
)
{
    for (int i = 0; i < Nsound; ++i)
    {
        // MATLAB: crit = (K_1_bar*M/(alpha*beta))
        Real crit = (K_1_bar[i] * M / (alpha * beta));

        Real Ireal = 0.0, Iimag = 0.0;

        if (K_2_bar[i] * K_2_bar[i] < (crit * crit))
        {
            // subcritical branch
            Real k_min_bar = std::sqrt(std::abs(
                mu_bar[i] * mu_bar[i] -
                (K_2_bar[i] * K_2_bar[i]) / (beta * beta)));

            Real B = K_1_bar[i] + M * mu_bar[i] + k_min_bar;

            // f1 (identical to Amiet)
            Real fr1, fi1;
            Radiation_integral1(B, C[i], fr1, fi1);

            // f2 (back-scattering in R&M)
            Real fr2, fi2;
            Radiation_integral2(
                B,              // B
                K_bar[i],       // K_bar
                k_min_bar,      // k_min_bar
                mu_bar[i],      // mu_bar
                S0,             // S0
                K_1_bar[i],     // K_1_bar
                alpha,          // alpha
                x,              // x
                M,              // M
                fr2, fi2);

            Ireal = fr1 + fr2;
            Iimag = fi1 + fi2;
        }
        else
        {
            // supercritical branch
            Real k_min_bar_prime;
            Real A1prime_r, A1prime_i;
            Real Aprime_r, Aprime_i;
            Real Thetaprime_r, Thetaprime_i;

            // compute wavenumbers etc.
            Wavenumbers_subcrit(
                M, mu_bar[i], K_1_bar[i], K_2_bar[i], K_bar[i], beta,
                k_min_bar_prime, A1prime_r, A1prime_i,
                Aprime_r, Aprime_i,
                Thetaprime_r, Thetaprime_i);

            // f1
            Real fr1, fi1;
            Radiation_integral1_subcrit(
                C[i],
                A1prime_r, A1prime_i,
                mu_bar[i],
                x,
                S0,
                k_min_bar_prime,
                fr1, fi1);

            // f2
            Real fr2, fi2;
            Radiation_integral2_subcrit(
                Aprime_r, Aprime_i,
                A1prime_r, A1prime_i,
                mu_bar[i],
                M,
                x,
                S0,
                k_min_bar_prime,
                alpha,
                K_bar[i],
                Thetaprime_r, Thetaprime_i,
                fr2, fi2);

            Ireal = fr1 + fr2;
            Iimag = fi1 + fi2;
        }

        // store |I|^2 directly
        I_abs2[i] = Ireal * Ireal + Iimag * Iimag;
    }
}



void TE_noise_outer(
    // Flow / geometry parameters formerly in 'inputs':
    Real M, Real U, Real x, Real y, Real z,
    Real b, Real Ky, Real c, Real span,

    // Fluid properties formerly in 'fluid':
    Real c0, Real rho, Real nu,

    // Array of frequencies:
    const Real omega[Nsound],
    Real (&WPS_lower)[Nsound], Real (&WPS_upper)[Nsound],

    // Output wall-pressure spectrum:
    Real (&farfieldSpectra)[Nsound],
    Real (&farfieldSpectra_up)[Nsound],
    Real (&farfieldSpectra_down)[Nsound]
)
{


    // Finding wavenumbers etc /////////////////////////////////
    Real beta = std::sqrt(1.0 - M*M);
    // distance to observer
    Real S0 = std::sqrt(x*x + beta*beta * (y*y + z*z));
    // convection velocity
    Real U_c = 0.7 * U; 
    Real alpha = U / U_c;

    // Outputs (arrays):
    Real C[Nsound];
    Real K_bar[Nsound];
    Real mu_bar[Nsound];
    Real K_1_bar[Nsound];
    Real K_2_bar[Nsound];
    Real k_bar[Nsound];
    for (int i=0; i<Nsound; ++i)
    {
        // acoustic wavenumber
        Real k = omega[i] / c0; // correct

        // aerodynamic wavenumber
        Real K = omega[i] / U; // correct

        // non-dimensional forms
        K_bar[i]    = K * b; //both correct
        k_bar[i]= k * b;

        // mu_bar
        mu_bar[i]   = K_bar[i] * M / (beta*beta); //correct

        // K_1_bar
        K_1_bar[i]  = alpha * K_bar[i]; //correct 

        // C (temporal variable)
        C[i] = K_1_bar[i] - mu_bar[i] * ((x/S0) - M); //correct

        // K_2_bar
        if (Ky == 0.0) {
            K_2_bar[i] = k_bar[i] * y / S0;
        } else {
            // If you really want the linspace(-Ky,Ky,100) branch,
            // you’d handle it separately. For now we mimic Ky==0 case.
            K_2_bar[i] = k_bar[i] * y / S0;
        }
    }


    // when Ky = 0 :
    Real I_abs2[Nsound];
    Radiation_integral_total(C,K_bar,mu_bar,S0,K_1_bar,K_2_bar,alpha,M,x,
                            beta,Ky,I_abs2);
    

    // spanwise correlation length :
    Real b_c = 1.47; // corcos constant
    Real K_2 = 0.0;
    
    Real l_y[Nsound];
    for (int i=0;i<Nsound;++i){
        
        Real top = omega[i] / (b_c*U_c);
        Real bot = (omega[i]*omega[i]) / (b_c*U_c * b_c*U_c);
        l_y[i] = top/bot;
    }


    // far field spectra (eq 18 in R&M) :
    for (int i=0;i<Nsound;++i){
        
        Real term1 = std::pow((omega[i]*c*z)/(c0*2.0*2.0*M_PI*S0*S0), 2.0);

        Real Spp_upper = term1*2.0*M_PI*span*I_abs2[i]*(WPS_upper[i]*l_y[i]/M_PI);
        Real Spp_lower = term1*2.0*M_PI*span*I_abs2[i]*(WPS_lower[i]*l_y[i]/M_PI);

        farfieldSpectra_up[i] = Spp_upper;
        farfieldSpectra_down[i] = Spp_lower;
        farfieldSpectra[i] = Spp_upper + Spp_lower ;
    }

}