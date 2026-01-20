#pragma once

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <codi.hpp>
#include "real_type.h"
#include "data_structs.h"

struct Isol;

template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
 
template<typename Active>
struct ImplicitBlockSolveData {
  using Real       = typename Active::Real;
  using Identifier = typename Active::Identifier;
  using Tape       = typename Active::Tape;
  using VectorAccess = codi::VectorAccessInterface<Real, Identifier>;

  int n; // Ncoords+1
  int m; // Npanels

  // Stored primal values (passive)
  Eigen::MatrixXd A;   // n x n
  Eigen::MatrixXd X;   // n x m   (solution, i.e. Bp_full)

  // Identifiers for inputs/outputs (needed to read/write adjoints)
  std::vector<Identifier> A_id;   // n*n
  std::vector<Identifier> B_id;   // n*m
  std::vector<Identifier> X_id;   // (Ncoords)*m  (your output storage size)

  ImplicitBlockSolveData(int n_, int m_)
    : n(n_), m(m_), A(n_, n_), X(n_, m_),
      A_id(n_*n_), B_id(n_*m_), X_id((n_-1)*m_) {}
};


// Reverse callback:
// Given adjoints of X (Bp) -> compute adjoints of A and B using implicit rule
template<typename Active>
static void implicit_block_solve_b(typename Active::Tape* /*tape*/, void* d,
                                   codi::VectorAccessInterface<typename Active::Real, typename Active::Identifier>* adj)
{
  using Real       = typename Active::Real;
  using Identifier = typename Active::Identifier;

  auto* data = static_cast<ImplicitBlockSolveData<Active>*>(d);

  const int n = data->n;
  const int m = data->m;

  const size_t maxDim = adj->getVectorSize();

  // Factorization of A^T (passive) - reuse per dimension
  Eigen::PartialPivLU<Eigen::MatrixXd> luAT(data->A.transpose());

  for (size_t dim = 0; dim < maxDim; ++dim) {

    // 1) Gather X_b (n x m) from tape adjoints.
    // Your output only stores first (n-1)=Ncoords rows; last row is treated as 0 adjoint.
    Eigen::MatrixXd X_b = Eigen::MatrixXd::Zero(n, m);

    for (int col = 0; col < m; ++col) {
      for (int row = 0; row < n-1; ++row) {
        const int k = colMajorIndex(row, col, n-1);
        Identifier id = data->X_id[k];
        X_b(row, col) = adj->getAdjoint(id, dim);
        adj->resetAdjoint(id, dim); // x_b = 0 in CoDiPack’s convention
      }
    }

    // 2) Solve adjoint system:  A^T * Lambda = X_b
    Eigen::MatrixXd Lambda = luAT.solve(X_b); // n x m

    // 3) Update adjoints:
    //    B_b += -Lambda   (because A X = -B)
    for (int col = 0; col < m; ++col) {
      for (int row = 0; row < n; ++row) {
        const int k = colMajorIndex(row, col, n);
        Identifier id = data->B_id[k];
        adj->updateAdjoint(id, dim, Real(-Lambda(row, col)));
      }
    }

    //    A_b += -(Lambda * X^T)
    //    (A is n x n, Lambda n x m, X n x m => Lambda*X^T is n x n)
    Eigen::MatrixXd dA = -(Lambda * data->X.transpose());

    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        const int k = colMajorIndex(i, j, n);
        Identifier id = data->A_id[k];
        adj->updateAdjoint(id, dim, Real(dA(i, j)));
      }
    }
  }
}

template<typename Active>
static void implicit_block_solve_delete(typename Active::Tape* /*tape*/, void* d)
{
  delete static_cast<ImplicitBlockSolveData<Active>*>(d);
}


template<typename Real>
void solve_sys_ue(Isol& isolc, const Real* RHS, Real* Bp)
{
  constexpr int n = Ncoords + 1;
  constexpr int m = Ncoords + Nwake - 2;  // nPanels

  using Active = Real;
  using Tape   = typename Active::Tape;

  Tape& tape = Active::getTape(); // this exists for CoDiPack ActiveType

  // Allocate external-function data (only used if tape is active)
  auto* data = new ImplicitBlockSolveData<Active>(n, m);

  // 1) Extract primal A (n x n) and store identifiers
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      const int k = colMajorIndex(i, j, n);
      const Active& a = isolc.infMatrix[k];
      data->A(i, j) = a.getValue();
      data->A_id[k] = a.getIdentifier();
    }
  }

  // 2) Extract primal B (n x m) and store identifiers
  Eigen::MatrixXd Bval(n, m);
  for (int col = 0; col < m; ++col) {
    for (int row = 0; row < n; ++row) {
      const int k = colMajorIndex(row, col, n);
      const Active& b = RHS[k];
      Bval(row, col) = b.getValue();
      data->B_id[k]  = b.getIdentifier();
    }
  }

  // 3) Passive solve: A * X = -B
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(data->A);
  data->X = lu.solve(-Bval); // n x m

  // 4) Write primal output values into Bp (Ncoords x m)
  for (int col = 0; col < m; ++col) {
    for (int row = 0; row < Ncoords; ++row) {
      const int outk = colMajorIndex(row, col, Ncoords);
      Bp[outk] = data->X(row, col); // primal assignment
    }
  }

  // 5) If tape inactive, we're done
  if (!tape.isActive()) {
    delete data;
    return;
  }

  // 6) Register outputs on tape + store their identifiers
  //    IMPORTANT: CoDiPack needs outputs registered as "external function outputs"
  for (int col = 0; col < m; ++col) {
    for (int row = 0; row < Ncoords; ++row) {
      const int outk = colMajorIndex(row, col, Ncoords);
      // registerExternalFunctionOutput mutates the variable's identifier bookkeeping
      (void)tape.registerExternalFunctionOutput(Bp[outk]);
      data->X_id[outk] = Bp[outk].getIdentifier();
    }
  }

  // 7) Push external function with reverse callback
  tape.pushExternalFunction(
    codi::ExternalFunction<Tape>::create(
      /*reverse*/ &implicit_block_solve_b<Active>,
      /*data*/    data,
      /*del*/     &implicit_block_solve_delete<Active>
      // forward/primal callbacks omitted => CoDiPack uses stored primals already
    )
  );
};