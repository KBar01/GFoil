#include <complex>
#include <cmath>
#include "Faddeeva.hh"

constexpr int N1 = 80;
constexpr int N2 = 86;
constexpr int N3 = 21;
constexpr int N4 = 13;


using cplx = std::complex<double>;

void fresnel(const cplx& z, cplx& S_out, cplx& C_out) {
    double w0 = std::abs(z);
    cplx zp = 0.5 * M_PI * z * z;
    cplx zp2 = zp * zp;

    int branch_i = 0;
    if (w0 > 0.0) branch_i++;
    if (w0 > 2.5) branch_i++;
    if (w0 >= 4.5) branch_i++;

    cplx S = 0.0, C = 0.0;

    if (branch_i == 0) {
        S = 0.0;
        C = 0.0;
    } 

    else if (branch_i == 1) {
        S = z * zp / 3.0;
        cplx cr = S;

        for (int k = 1; k <= N1; ++k) {
            cr *= -0.5 * (4.0 * k - 1.0) / (k * (2.0 * k + 1.0) * (4.0 * k + 3.0)) * zp2;
            S += cr;
        }

        C = z;
        cr = C;

        for (int k = 1; k <= N1; ++k) {
            cr *= -0.5 * (4.0 * k - 3.0) / (k * (2.0 * k - 1.0) * (4.0 * k + 1.0)) * zp2;
            C += cr;
        }
    } 
    
    else if (branch_i == 2) {
        cplx cf0 = 1e-100;
        cplx cf1 = 0.0;
        cplx cf;
        S = 0.0;
        C = 0.0;

        for (int i = N2 - 1; i >= 0; --i) {
            int k = i;
            cf = (2.0 * k + 3.0) * cf0 / zp - cf1;
            if (k % 2 == 0)
                C += cf;
            else
                S += cf;
            cf1 = cf0;
            cf0 = cf;
        }

        cplx coeff = 2.0 / (M_PI * z) * std::sin(zp) / cf0;
        S *= coeff;
        C *= coeff;
    } 
    
    else {
        double x = z.real();
        double y = z.imag();

        cplx DS, DC;
        if (y < x) {
            DS = (y < -x) ? cplx(0.0, 0.5) : cplx(0.5, 0.0);
            DC = (y < -x) ? cplx(0.0, -0.5) : cplx(0.5, 0.0);
        } else {
            DS = (y < -x) ? cplx(0.0, -0.5) : cplx(-0.5, 0.0);
            DC = (y < -x) ? cplx(0.0, 0.5)  : cplx(-0.5, 0.0);
        }

        cplx cf = 1.0, cr = 1.0;
        for (int k = 1; k <= N3; ++k) {
            cr *= -0.25 * (4.0 * k - 1.0) * (4.0 * k - 3.0) / zp2;
            cf += cr;
        }

        cplx cg_S = 1.0, cr_S = 1.0;
        for (int k = 1; k <= N4; ++k) {
            cr_S *= -0.25 * (4.0 * k + 1.0) * (4.0 * k - 1.0) / zp2;
            cg_S += cr_S;
        }
        cg_S /= (M_PI * z * z);

        cplx cg_C = 1.0 / (M_PI * z * z), cr_C = cg_C;
        for (int k = 1; k <= N4; ++k) {
            cr_C *= -0.25 * (4.0 * k + 1.0) * (4.0 * k - 1.0) / zp2;
            cg_C += cr_C;
        }

        S = DS - (cf * std::cos(zp) + cg_S * std::sin(zp)) / (M_PI * z);
        C = DC + (cf * std::sin(zp) - cg_C * std::cos(zp)) / (M_PI * z);
    }

    S_out = S;
    C_out = C;
}


void error_func(const cplx& z, cplx& z_out){
    
    z_out = Faddeeva::erf(z) ;

}

