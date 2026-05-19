#pragma once

#include <vector>
#include <string>
#include <array>
#include "real_type.hpp"

// Nine shared structs (field-identical with fwd data_structs.h) from the
// shared template header; the template aliases below keep the same names.
#include "data_structs_shared.hpp"

template<typename Real> using Geom  = Geom_t<Real>;
template<typename Real> using TE    = TE_t<Real>;
template<typename Real> using Foil  = Foil_t<Real>;
template<typename Real> using Wake  = Wake_t<Real>;
template<typename Real> using Oper  = Oper_t<Real>;
template<typename Real> using Vsol  = Vsol_t<Real>;
template<typename Real> using Post  = Post_t<Real>;
template<typename Real> using Param = Param_t<Real>;
// Trans alias removed: forced transition deleted; natural e^9 only.

//-------------------------------------------------------------------------------
// Inviscid Solution Struct — split into constant (c) and variable (v) parts
template<typename Real>
struct Isolc {
    Real infMatrix[(Ncoords + 1) * (Ncoords + 1)] = {0};
    Real gammasRef[2 * Ncoords] = {0};
    Real gammas[Ncoords] = {0};
    Real uewi[Nwake]     = {0};
    Real uewiref[2 * Nwake] = {0};
};

template<typename Real>
struct Isolv {
    Real stagArcLocation  = 0;
    Real stagXLocation[2] = {0};
    Real sstag_ue[2]      = {0};
    int  stagIndex[2]     = {0};
    int  edgeVelSign[Ncoords] = {0};
    Real distFromStag[Ncoords + Nwake] = {0};

    Isolv() {
        for (int i = 0; i < Ncoords; ++i)
            edgeVelSign[i] = -1;
    }
};

//-------------------------------------------------------------------------------
// Global solve struct (AD: no sparse Jacobian storage — CoDi handles derivatives)
template<typename Real>
struct Glob {
    Real U[4*(Ncoords+Nwake)] = {0};
    Real R[4*(Ncoords+Nwake)] = {0.0};
};
