#pragma once

#include <codi.hpp>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>
#include "real_type.hpp"
#include "data_structs.hpp"
#include "get_funcs.hpp"
#include "panel_funcs.hpp"


template<typename Real> struct Isolc;
template<typename Real> struct Isolv;
template<typename Real> struct Vsol;
template<typename Real> struct Foil;
template<typename Real> struct Param;
template<typename Real> struct Post;
template<typename Real> struct Oper;
template<typename Real> struct Geom;
template<typename Real> struct Wake;
template<typename Real> struct Glob;
template<typename Real> struct Trans;


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

template<typename Real>
void solve_sys_inv(Isolc<Real> &isol, const Real* RHS)
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


template<typename Real>
void build_gamma_codi(Isolc<Real> &isol, const Foil<Real>& foil, const Oper<Real>& op) {
    // Build and solve the inviscid linear system for alpha=0,90,input
    // INPUT
    //   isol: Invisid solution structure
    //   foil: Aerofoil data  structure
    //   alpha: angle of attack (Radians)
    // OUTPUT
    //   isol.gamref: 0,90deg vorticity distributions at each node (Nx2)
    //   isol.gam: gamma for the particular input angle, alpha
    //   isol.infMatrix: aerodynamic influence coefficient matrix, filled in

    // Initialize matrices A and rhs
    Real rhs[2*(Ncoords + 1)] ; 

    // Initialise the panelinfo struct and influence coeffs used inside loop
    PanelInfo<Real> panelInfo;
    Real aij = 1.0;
    Real bij = 1.0;
    Real a   = 1.0;
    Real a_vortex = 1.0;
    Real b_vortex = 1.0;

    // Build influence matrix and rhs
    for (int i = 0; i < Ncoords; i++) {  // Loop over nodes

        for (int j = 0; j < Ncoords-1; j++) {  // Loop over panels
            // panel_linvortex_stream should be implemented to return aij and bij
            panel_linvortex_stream(foil.x[colMajorIndex(0,j,2)], foil.x[colMajorIndex(1,j,2)],
            foil.x[colMajorIndex(0,j+1,2)], foil.x[colMajorIndex(1,j+1,2)],
            foil.x[colMajorIndex(0,i,2)], foil.x[colMajorIndex(1,i,2)],
            panelInfo, aij, bij);
            
            isol.infMatrix[colMajorIndex(i,j,Ncoords+1)]   += aij;
            isol.infMatrix[colMajorIndex(i,j+1,Ncoords+1)] += bij;
        }

        isol.infMatrix[colMajorIndex(i,Ncoords,Ncoords+1)] = -1.0;  // Last unknown = streamfunction value on surface
        
        // Right-hand sides
        rhs[colMajorIndex(i,0,Ncoords+1)] = -foil.x[colMajorIndex(1,i,2)];
        rhs[colMajorIndex(i,1,Ncoords+1)] = foil.x[colMajorIndex(0,i,2)];

        // TE source influence
        panel_constsource_stream(foil.x[colMajorIndex(0,Ncoords-1,2)], foil.x[colMajorIndex(1,Ncoords-1,2)],
                              foil.x[colMajorIndex(0,0,2)], foil.x[colMajorIndex(1,0,2)],
                              foil.x[colMajorIndex(0,i,2)], foil.x[colMajorIndex(1,i,2)],
                              panelInfo, a);
        
        isol.infMatrix[colMajorIndex(i,0,Ncoords+1)] += -a*(0.5 * foil.te.tcp);
        isol.infMatrix[colMajorIndex(i,Ncoords-1,Ncoords+1)] += a*(0.5 * foil.te.tcp);

        // TE vortex panel 
        panel_linvortex_stream(foil.x[colMajorIndex(0,Ncoords-1,2)], foil.x[colMajorIndex(1,Ncoords-1,2)],
        foil.x[colMajorIndex(0,0,2)], foil.x[colMajorIndex(1,0,2)],
        foil.x[colMajorIndex(0,i,2)], foil.x[colMajorIndex(1,i,2)],
        panelInfo, a_vortex, b_vortex);
        
        isol.infMatrix[colMajorIndex(i,0,Ncoords+1)] += -(a_vortex + b_vortex)*(-0.5 * foil.te.tdp);
        isol.infMatrix[colMajorIndex(i,Ncoords-1,Ncoords+1)] += (a_vortex + b_vortex)*(-0.5 * foil.te.tdp);
    }

    // Kutta condition
    isol.infMatrix[colMajorIndex(Ncoords, 0, Ncoords+1)] = 1;
    isol.infMatrix[colMajorIndex(Ncoords, Ncoords-1, Ncoords+1)] = 1;
    //solve_sys(isol,rhs);   // Solves the linear sys using Eigen, and puts gammas in isol struct

    
    solve_sys_inv(isol,rhs);


    for (int i = 0; i < Ncoords; ++i) {
        isol.gammas[i] = isol.gammasRef[colMajorIndex(i,0,Ncoords)]*std::cos(op.alpha) + isol.gammasRef[colMajorIndex(i,1,Ncoords)]*std::sin(op.alpha);
    }
};
