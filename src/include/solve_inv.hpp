#pragma once

// Shared inviscid-solve header — used by both GFoil_fwd_codi and GFoil_AD.
//
// solve_sys_inv : CoDiPack-aware linear solve for the AIC system.
//   typename Real  — CoDiPack active type (deduced from RHS array)
//   typename IsolT — duck-typed: accepts Isol (fwd) or Isolc<Real> (AD)
//                    both provide infMatrix[] and gammasRef[].
//
// build_gamma_codi : assemble AIC matrix and call solve_sys_inv.
//   typename IsolT, FoilT, OperT — duck-typed for the same reason.
//
// Callers must include the appropriate real_type and data_structs headers
// before including this header.  panel_funcs.hpp is included here because
// the panel geometry functions are needed inside build_gamma_codi.

#include <codi.hpp>
#include <Eigen/Dense>
#include <cmath>
#include "panel_funcs.hpp"

// ── CoDiPack external-function data and callbacks (identical for fwd and AD) ──

template<typename Active>
struct ImplicitInvSolveData {
  using Real       = typename Active::Real;
  using Identifier = typename Active::Identifier;
  using Tape       = typename Active::Tape;

  int n; // Ncoords+1

  Eigen::MatrixXd A;   // n x n
  Eigen::MatrixXd X;   // n x 2 (solutions)

  std::vector<Identifier> A_id; // n*n
  std::vector<Identifier> B_id; // n*2
  std::vector<Identifier> X_id; // (Ncoords)*2

  ImplicitInvSolveData(int n_)
    : n(n_), A(n_,n_), X(n_,2),
      A_id(n_*n_), B_id(n_*2), X_id((n_-1)*2) {}
};

template<typename Active>
static void implicit_inv_solve_b(typename Active::Tape*,
                                 void* d,
                                 codi::VectorAccessInterface<
                                   typename Active::Real,
                                   typename Active::Identifier>* adj)
{
  using Real = typename Active::Real;

  auto* data = static_cast<ImplicitInvSolveData<Active>*>(d);
  const int n = data->n;

  Eigen::PartialPivLU<Eigen::MatrixXd> luAT(data->A.transpose());

  const size_t maxDim = adj->getVectorSize();

  for (size_t dim = 0; dim < maxDim; ++dim) {

    // 1) Collect X_b (n x 2)
    Eigen::MatrixXd X_b = Eigen::MatrixXd::Zero(n, 2);

    for (int col = 0; col < 2; ++col) {
      for (int row = 0; row < n-1; ++row) {
        int k = colMajorIndex(row, col, n-1);
        auto id = data->X_id[k];
        X_b(row, col) = adj->getAdjoint(id, dim);
        adj->resetAdjoint(id, dim);
      }
    }

    // 2) Adjoint solve: A^T Λ = X_b
    Eigen::MatrixXd Lambda = luAT.solve(X_b); // n x 2

    // 3) B adjoint update: B_b += Lambda
    for (int col = 0; col < 2; ++col) {
      for (int row = 0; row < n; ++row) {
        int k = colMajorIndex(row, col, n);
        adj->updateAdjoint(data->B_id[k], dim,
                           Real(Lambda(row,col)));
      }
    }

    // 4) A adjoint update: A_b += -(Lambda X^T)
    Eigen::MatrixXd dA = -(Lambda * data->X.transpose());
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        int k = colMajorIndex(i,j,n);
        adj->updateAdjoint(data->A_id[k], dim,
                           Real(dA(i,j)));
      }
  }
}

template<typename Active>
static void implicit_inv_solve_delete(typename Active::Tape*, void* d)
{
  delete static_cast<ImplicitInvSolveData<Active>*>(d);
}

// ── solve_sys_inv ─────────────────────────────────────────────────────────────
// IsolT accepts Isol (forward) or Isolc<Real> (AD) — both expose infMatrix[]
// and gammasRef[].  Real is deduced from the RHS pointer type.

template<typename Real, typename IsolT>
void solve_sys_inv(IsolT& isol, const Real* RHS)
{
  constexpr int n = Ncoords + 1;
  using Active = Real;
  using Tape   = typename Active::Tape;

  Tape& tape = Active::getTape();

  auto* data = new ImplicitInvSolveData<Active>(n);

  // --- Extract A ---
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      int k = colMajorIndex(i,j,n);
      const Active& a = isol.infMatrix[k];
      data->A(i,j)   = a.getValue();
      data->A_id[k]  = a.getIdentifier();
    }

  // --- Extract RHS into B (n x 2) ---
  Eigen::MatrixXd Bval(n,2);
  for (int i = 0; i < n; ++i) {
    const Active& b1 = RHS[i];
    const Active& b2 = RHS[n + i];

    Bval(i,0) = b1.getValue();
    Bval(i,1) = b2.getValue();

    data->B_id[colMajorIndex(i,0,n)] = b1.getIdentifier();
    data->B_id[colMajorIndex(i,1,n)] = b2.getIdentifier();
  }

  // --- Passive solve ---
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(data->A);
  data->X = lu.solve(Bval); // n x 2

  // --- Write primal outputs ---
  for (int i = 0; i < Ncoords; ++i) {
    isol.gammasRef[i]           = data->X(i,0);
    isol.gammasRef[Ncoords + i] = data->X(i,1);
  }

  if (!tape.isActive()) {
    delete data;
    return;
  }

  // --- Register outputs ---
  for (int col = 0; col < 2; ++col)
    for (int row = 0; row < Ncoords; ++row) {
      int k = colMajorIndex(row, col, Ncoords);
      Real& out = isol.gammasRef[col*Ncoords + row];
      tape.registerExternalFunctionOutput(out);
      data->X_id[k] = out.getIdentifier();
    }

  // --- Push external function ---
  tape.pushExternalFunction(
    codi::ExternalFunction<Tape>::create(
      &implicit_inv_solve_b<Active>,
      data,
      &implicit_inv_solve_delete<Active>
    )
  );
}

// ── build_gamma_codi ──────────────────────────────────────────────────────────
// Assembles the AIC matrix, builds the RHS, and calls solve_sys_inv.
// IsolT, FoilT, OperT are duck-typed: works with the non-template
// (Isol, Foil, Oper) forward structs and the template (Isolc<Real>, Foil<Real>,
// Oper<Real>) AD structs provided they expose the same member names.

template<typename Real, typename IsolT, typename FoilT, typename OperT>
void build_gamma_codi(IsolT& isol, const FoilT& foil, const OperT& op)
{
    Real rhs[2*(Ncoords + 1)];

    // Precompute panel-fixed geometry (t, n, d) once per panel.
    // Indices 0..Ncoords-2: body panels (j → j+1).
    // Index Ncoords-1: TE panel (Ncoords-1 → 0).
    PanelGeom<Real> panelGeoms[Ncoords];
    for (int j = 0; j < Ncoords - 1; ++j) {
        precompute_panel_geom<Real>(
            foil.x[colMajorIndex(0,j,2)],   foil.x[colMajorIndex(1,j,2)],
            foil.x[colMajorIndex(0,j+1,2)], foil.x[colMajorIndex(1,j+1,2)],
            panelGeoms[j]);
    }
    precompute_panel_geom<Real>(
        foil.x[colMajorIndex(0,Ncoords-1,2)], foil.x[colMajorIndex(1,Ncoords-1,2)],
        foil.x[colMajorIndex(0,0,2)],          foil.x[colMajorIndex(1,0,2)],
        panelGeoms[Ncoords-1]);

    PanelInfo<Real> panelInfo;
    Real aij = 1.0;
    Real bij = 1.0;
    Real a   = 1.0;
    Real a_vortex = 1.0;
    Real b_vortex = 1.0;

    for (int i = 0; i < Ncoords; i++) {

        for (int j = 0; j < Ncoords-1; j++) {
            panel_linvortex_stream(
                panelGeoms[j],
                foil.x[colMajorIndex(0,j,2)], foil.x[colMajorIndex(1,j,2)],
                foil.x[colMajorIndex(0,i,2)], foil.x[colMajorIndex(1,i,2)],
                panelInfo, aij, bij);

            isol.infMatrix[colMajorIndex(i,j,Ncoords+1)]   += aij;
            isol.infMatrix[colMajorIndex(i,j+1,Ncoords+1)] += bij;
        }

        isol.infMatrix[colMajorIndex(i,Ncoords,Ncoords+1)] = -1.0;

        rhs[colMajorIndex(i,0,Ncoords+1)] = -foil.x[colMajorIndex(1,i,2)];
        rhs[colMajorIndex(i,1,Ncoords+1)] =  foil.x[colMajorIndex(0,i,2)];

        panel_constsource_stream(
            panelGeoms[Ncoords-1],
            foil.x[colMajorIndex(0,Ncoords-1,2)], foil.x[colMajorIndex(1,Ncoords-1,2)],
            foil.x[colMajorIndex(0,i,2)],          foil.x[colMajorIndex(1,i,2)],
            panelInfo, a);

        isol.infMatrix[colMajorIndex(i,0,Ncoords+1)]         += -a*(0.5 * foil.te.tcp);
        isol.infMatrix[colMajorIndex(i,Ncoords-1,Ncoords+1)] +=  a*(0.5 * foil.te.tcp);

        panel_linvortex_stream(
            panelGeoms[Ncoords-1],
            foil.x[colMajorIndex(0,Ncoords-1,2)], foil.x[colMajorIndex(1,Ncoords-1,2)],
            foil.x[colMajorIndex(0,i,2)],          foil.x[colMajorIndex(1,i,2)],
            panelInfo, a_vortex, b_vortex);

        isol.infMatrix[colMajorIndex(i,0,Ncoords+1)]         += -(a_vortex+b_vortex)*(-0.5*foil.te.tdp);
        isol.infMatrix[colMajorIndex(i,Ncoords-1,Ncoords+1)] +=  (a_vortex+b_vortex)*(-0.5*foil.te.tdp);
    }

    // Kutta condition
    isol.infMatrix[colMajorIndex(Ncoords, 0,         Ncoords+1)] = 1;
    isol.infMatrix[colMajorIndex(Ncoords, Ncoords-1, Ncoords+1)] = 1;

    solve_sys_inv(isol, rhs);

    for (int i = 0; i < Ncoords; ++i) {
        isol.gammas[i] = isol.gammasRef[colMajorIndex(i,0,Ncoords)] * std::cos(op.alpha)
                       + isol.gammasRef[colMajorIndex(i,1,Ncoords)] * std::sin(op.alpha);
    }
}


// Fwd non-template wrapper — only when data_structs.h has been included.
#ifdef AIRFOIL_STRUCTS_H
inline void build_gamma_codi(Isol& isol, const Foil& foil, const Oper& op) {
    build_gamma_codi<Real>(isol, foil, op);
}
#endif // AIRFOIL_STRUCTS_H
