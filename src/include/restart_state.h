#pragma once
#include <string>
#include <vector>

struct RestartState {
    std::vector<double> states;
    std::vector<int>    turb;
    std::vector<int>    stag;
    std::vector<double> RVvals;
    std::vector<int>    RVrows;
    std::vector<int>    RVcols;
    int                 RVnz = 0;
};

struct ForwardResult {
    bool   converged = false;
    double CL = 0.0, CD = 0.0, CM = 0.0, OASPL = 0.0;
    // Empty when converged; one of "transition_front_oscillation", "diverged",
    // "no_convergence" when not converged.
    std::string failure_mode = "";

    // ── verbose output (only populated when verbose=true) ────────────────────
    std::vector<double> innerFoilX;   // length Ncoords
    std::vector<double> innerFoilY;   // length Ncoords

    std::vector<double> Cp;           // pressure coefficient
    std::vector<double> delta_star;   // displacement thickness [m]
    std::vector<double> theta;        // momentum thickness [m]
    std::vector<double> tau_wall;     // wall shear stress [Pa]
    std::vector<double> tau_max;      // max shear stress [Pa]; 0 if laminar
    std::vector<double> Ue;           // BL edge velocity [m/s]
    std::vector<double> dpdx;         // streamwise pressure gradient [Pa/m]
    std::vector<bool>   is_turb;      // turbulence flag per node

    double topTransX = 0.0;
    double botTransX = 0.0;

    // [theta, delta*, tau_max, Ue, dpdx, tau_wall, delta99]
    std::vector<double> BL_top;       // length 7
    std::vector<double> BL_bot;       // length 7

    std::vector<double> freq_Hz;      // length Nsound
    std::vector<double> WPS_upper;    // wall-pressure PSD upper [Pa^2/Hz]
    std::vector<double> WPS_lower;    // wall-pressure PSD lower [Pa^2/Hz]
    int nObs = 0;
    std::vector<double> FF_spectra;   // flat row-major (nObs, Nsound)

    int newton_iterations = 0;        // converging Newton iteration (param.niglob if not converged)
};
