
#pragma once

#include <codi.hpp>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>
#include "real_type.h"
#include "data_structs.h"


struct Isol ;

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
void solve_sys_inv(Isol&isol, const Real* RHS)
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
};