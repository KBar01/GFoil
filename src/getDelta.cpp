#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "real_type.h"
#include "main_func.h"
// Spalding law of the wall constants
const Real kappa = 0.41;
const Real A = 0.1108;

// Numerically invert Spalding's law for uplus given yplus
Real solve_spalding(Real yplus, Real tol = 1e-6, int max_iter = 30) {
    Real uplus = yplus < 11.0 ? yplus : std::log(yplus) / kappa; // Initial guess
    for (int i = 0; i < max_iter; ++i) {
        Real expkU = std::exp(kappa * uplus);
        Real F = yplus - uplus - A * (expkU - 1.0 - kappa * uplus
                  - 0.5 * std::pow(kappa * uplus, 2) - (1.0 / 6.0) * std::pow(kappa * uplus, 3));
        // Numerical derivative
        Real h = 1e-6;
        Real expkU_h = std::exp(kappa * (uplus + h));
        Real Fh = yplus - (uplus + h) - A * (expkU_h - 1.0 - kappa * (uplus + h)
                  - 0.5 * std::pow(kappa * (uplus + h), 2) - (1.0 / 6.0) * std::pow(kappa * (uplus + h), 3));
        Real dF = (Fh - F) / h;
        if (std::abs(dF) < 1e-10) dF = -1.0; // Prevent divide by near-zero
        Real du = -F / dF;
        uplus += du;
        if (std::abs(du) < tol) break;
    }
    return std::max(0.0, uplus); // Should not go negative
}

// Generate turbulent BL profile, outputting y, u/Ue, and the calculated delta
void turbulent_BL_profile_XFOIL(
    Real theta,        // Momentum thickness [m]
    Real Ue,           // Edge velocity [m/s]
    Real utau,         // Friction velocity [m/s]
    Real nu,           // Kinematic viscosity [m^2/s]
    Real& delta                    // Output: BL thickness (where u/Ue = 0.99)
) {


    constexpr int N = 100;
    const Real Nreal = 100.0;
    // Typical max y to search for edge (should be much larger than theta)
    Real ymax = 12.0 * theta;

    Real y[N];
    Real u_over_Ue[N];

    // Power law exponent
    const Real power = 1.0 / 7.0;

    // We will record where u/Ue first exceeds 0.99
    delta = -1.0;
    for (int i = 0; i < N; ++i) {
        y[i] = ymax * i / (Nreal - 1);
        Real yplus = y[i] * utau / nu;
        Real u_y = 0.0;
        if (yplus < 100) {
            // Inner region: Spalding law
            Real uplus = solve_spalding(yplus);
            u_y = uplus * utau;
        } else {
            // Outer region: power law, matched at yplus=100
            Real y_match = 100.0 * nu / utau;
            Real uplus_match = solve_spalding(100.0);
            Real u_match = uplus_match * utau;
            Real u_outer = Ue * std::pow(y[i] / ymax, power);
            Real u_outer_match = Ue * std::pow(y_match / ymax, power);
            u_y = u_outer + (u_match - u_outer_match);
            u_y = std::min(u_y, Ue);
        }
        u_over_Ue[i] = u_y / Ue;
        u_over_Ue[i] = std::max(0.0, std::min(1.0, u_over_Ue[i]));

        // Find delta (where u/Ue >= 0.99), if not already set
        if (delta < 0.0 && u_over_Ue[i] >= 0.99) {
            delta = y[i];
        }
    }
    // If profile never reaches 0.99, set delta to max y
    if (delta < 0.0) delta = ymax;
}