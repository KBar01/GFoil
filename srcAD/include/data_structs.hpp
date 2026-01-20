#pragma once

#include <vector>
#include <string>
#include <array>
#include "real_type.hpp"

//-------------------------------------------------------------------------------
// Geometry Struct
template<typename Real>
struct Geom {

    // Constants in coupled solve (X)
    Real chord = 1.0;                  // Chord length
    Real wakelen = 1.0;                // Wake extent length (in chords)
    int npoint = 1;                       // Number of geometry representation points
    Real xref[2] = {0.25, 0};     // Moment reference point
};

//-------------------------------------------------------------------------------
// TE Struct
template<typename Real>
struct TE {

    // Constants in coupled solve (X)
    Real t[2];            // Bisector vector (normalized)
    Real hTE = 0;             // Trailing edge gap
    Real dtdx = 0;            // Thickness slope
    Real tcp = 0;             // |t cross p|
    Real tdp = 0;             // t dot p

    // Constructor for TE, which computes TE properties
    TE(const Real* coords) {
        compute_TE_info(coords);
    }

    // Function to compute TE properties
    void compute_TE_info(const Real* coords) {
        // Compute lower tangent vector
        Real lowerTang[2] = {coords[0] - coords[2], coords[1] - coords[3]};
        Real normLowerTang = norm2(lowerTang);
        lowerTang[0] /= normLowerTang;
        lowerTang[1] /= normLowerTang;

        // Compute upper tangent vector
        Real upperTang[2] = {coords[2 * (Ncoords - 1)] - coords[2 * (Ncoords - 2)], coords[2 * (Ncoords - 1) + 1] - coords[2 * (Ncoords - 2) + 1]};
        Real normUpperTang = norm2(upperTang);
        upperTang[0] /= normUpperTang;
        upperTang[1] /= normUpperTang;

        // Compute bisector vector
        t[0] = 0.5 * (lowerTang[0] + upperTang[0]);
        t[1] = 0.5 * (lowerTang[1] + upperTang[1]);
        Real t_norm = norm2(t);
        t[0] /= t_norm;
        t[1] /= t_norm;

        // Compute lower-to-upper connector vector
        Real TEVector[2] = {coords[2 * (Ncoords - 1)] - coords[0], coords[2 * (Ncoords - 1) + 1] - coords[1]};

        // Compute TE gap
        hTE = -TEVector[0] * t[1] + TEVector[1] * t[0];

        // Compute thickness slope
        dtdx = lowerTang[0] * upperTang[1] - upperTang[0] * lowerTang[1];

        // Compute unit vector along TE panel
        Real norm_s = norm2(TEVector);
        Real p[2] = {TEVector[0] / norm_s, TEVector[1] / norm_s};

        // Compute tcp and tdp
        tcp = std::abs(t[0] * p[1] - t[1] * p[0]);  // |t cross p|
        tdp = t[0] * p[0] + t[1] * p[1];  // t dot p
    }
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Foil {

    // Constants in coupled solve (X)
    Real x[2*Ncoords] = {0};   // Flattened 2D array (size: 2*N)
    Real s[Ncoords] = {0};       // Arc length values at nodes
    TE<Real> te;                 // Trailing edge properties

    Foil(const Real* coordinates)
        : te(coordinates) {
        // Copy coordinates into internal array
        for (int i = 0; i < 2*Ncoords; i++) {
            x[i] = coordinates[i];
        }
        compute_arclengths();
    }

    void compute_arclengths() {
        s[0] = 0.0;  // First point has zero arclength
        for (int i = 1; i < Ncoords; ++i) {
            Real dx = x[2*i] - x[2*(i-1)];
            Real dy = x[2* i + 1] - x[2 * (i-1) + 1];
            s[i] = s[i-1] + std::sqrt(dx*dx + dy*dy);
        }
    }
};

//-------------------------------------------------------------------------------
template<typename Real>
struct Wake {   
    
    // Constants in coupled solve (X)
    Real x[Nwake*2] = {0};       // Node coordinates (2xN) flattened
    Real s[Nwake] = {0};         // Arc length values at nodes
    Real t[Nwake*2] = {0};       // Tangents at nodes (dx/ds, dy/ds)
};


//-------------------------------------------------------------------------------
template<typename Real>
struct Oper {
    // Constants in coupled solve (X)

    Real Vinf = 1.0;          // Velocity magnitude
    Real alpha = 0.0;         // Angle of attack (degrees)
    Real rho = 1.225;           // Density
    Real Re = 1e5;            // Reynolds number
    Real Ma = 0.0;            // Mach number

    // Constructor for initialization
    Oper(Real alpha_in, Real Re_in, Real Ma_in) 
        : alpha(alpha_in), Re(Re_in), Ma(Ma_in) {}
};

//-------------------------------------------------------------------------------
// Inviscid Solution Struct
template<typename Real>
struct Isolc {
    
    // Constants in coupled solve (X)
    Real infMatrix[(Ncoords + 1) * (Ncoords + 1)] = {0}; // Aero influence coefficient matrix (flattened)
    Real gammasRef[2 * (Ncoords)] = {0};               // Vortex strengths at airfoil nodes (0,90-deg alpha)
    Real gammas[Ncoords] = {0};                      // Vortex strengths at airfoil nodes (current alpha)
    Real uewi[Nwake] = {0};                       // Inviscid edge velocity in wake
    Real uewiref[2 * Nwake] = {0};                // 0,90-deg alpha inviscid ue solutions on wake
};

template<typename Real>
struct Isolv {

    // Internal Variables in coupled Solve
    Real stagArcLocation = 0;                             // Stagnation point location
    Real stagXLocation[2] = {0};                      // Sensitivity of sstag w.r.t. adjacent gammas
    Real sstag_ue[2]= {0};                       // Sensitivity of sstag w.r.t. adjacent ue values
    int stagIndex[2] = {0};                             // Node indices before/after stagnation
    int edgeVelSign[Ncoords] = {0};                       // +1/-1 for upper/lower surface nodes
    Real distFromStag[Ncoords + Nwake] = {0};               // Distance from stagnation at all points
    

    // Constructor to initialize arrays
    Isolv() {
        // Initialize edgeVelSign to all 1
        for (int i = 0; i < Ncoords; ++i) {
            edgeVelSign[i] = -1;
        }
    }

};




//-------------------------------------------------------------------------------
// Viscous Solution Struct
template<typename Real>
struct Vsol {
    // Constants in coupled solve (X)
    Real wgap[Nwake] = {0.0};       // Wake gap over wake points
    Real ue_sigma[(Ncoords+Nwake)*(Ncoords+Nwake-2)] = {0};// d(ue)/d(source) matrix

    // Internal Variables in coupled Solve
    Real ue_m[(Ncoords+Nwake)*(Ncoords+Nwake)] = {0.0};// Linearization of ue w.r.t. mass (all nodes)
    
    // Internal non-tracked Variables in coupled solve
    bool turb[Ncoords+Nwake] = {false};// Flag (1 = turbulent, 0 = laminar)
    std::vector<std::vector<int>> Is; // Holds lower, upper, wake node indices
};

//-------------------------------------------------------------------------------
// Global Parameters Struct
template<typename Real>
struct Glob {

    Real U[4*(Ncoords+Nwake)] = {0};  // Primary states (th, ds, sa, ue)
    Real R[4*(Ncoords+Nwake)] = {0.0};                 // Residuals
};



//-------------------------------------------------------------------------------
// Post-Processing Struct
template<typename Real>
struct Post {

    Real cp[Ncoords+Nwake]= {0};       // Cp distribution
    //Real cpi[Ncoords] = {0};      // Inviscid Cp distribution
    Real cl = 0.0;              // Lift coefficient
    Real cd = 0.0;              // Total drag coefficient
    Real cm = 0.0;
};

//-------------------------------------------------------------------------------
// Parameters Struct
template<typename Real>
struct Param {

    // Constants in coupled solve (X)
    Real rtol = 1e-10;           // Residual tolerance
    int niglob = 50;               // Max global iterations
    // Viscous parameters
    Real ncrit = 9.0;
    Real Cuq = 1.0;
    Real Dlr = 0.9;
    Real SlagK = 5.6;
    // Initial Ctau after transition
    Real CtauC = 1.8;
    Real CtauE = 3.3;
    // G-beta locus parameters
    Real GA = 6.7, GB = 0.75, GC = 18.0;
    // Operating conditions
    Real Minf = 0.0, Vinf = 0.0, muinf = 0.0, mu0 = 0.0;
    Real rho0 = 1.0, H0 = 0.0, Tsrat = 0.35, gam = 1.4;
    Real KTb = 1.0, KTl = 0.0, cps = 0.0;
    int breakLoop = 100;
};

template<typename Real>
struct Trans {
    // Constants in coupled solve (X)
    int isForced[2] = {0,0}; // tracking if transition is forced or natural
    int transNode[2] = {1,198}; // start node of trans of bottom/top surface
    Real transPos[2] = {0.5,0.5}; // x position of trans of bottom/top surface
};
