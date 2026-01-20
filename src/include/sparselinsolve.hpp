#pragma once
#ifdef USE_CODIPACK
#include <codi.hpp>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <vector>
#include "real_type.h"
#include "data_structs.h"

// You already have these
inline int colMajorIndex(int i, int j, int n) { return i + j * n; };

// -----------------------------
// External-function data
// -----------------------------
template<typename Active>
struct ImplicitSparseSolveData {
  using Real       = typename Active::Real;
  using Identifier = typename Active::Identifier;
  using Tape       = typename Active::Tape;

  int n;      // system size
  int nnz;    // number of stored entries (glob.R_V_latest)

  // Passive data needed in reverse
  Eigen::SparseMatrix<double> A; // n x n
  Eigen::VectorXd x;             // solution of A x = b (NOT including your minus sign)

  // Store structure for each entry k
  std::vector<int> row;
  std::vector<int> col;

  // Identifiers (what we update in reverse)
  std::vector<Identifier> Aval_id; // size nnz: ids of glob.R_V_vals[k]
  std::vector<Identifier> b_id;    // size n:   ids of glob.R[i]
  std::vector<Identifier> out_id;  // size n:   ids of glob.dU[i] outputs

  ImplicitSparseSolveData(int n_, int nnz_)
    : n(n_), nnz(nnz_), A(n_, n_), x(n_),
      row(nnz_), col(nnz_),
      Aval_id(nnz_), b_id(n_), out_id(n_) {}
};

// -----------------------------
// Reverse callback
// -----------------------------
// Primal relationship implemented by the wrapper:
//    A x = b
//    dU  = -x
//
// Reverse:
//    x_b = -(dU_b)
//    Solve A^T lambda = x_b
//    b_b += lambda
//    Aval_b[k] += -(lambda[row_k] * x[col_k])
//
template<typename Active>
static void implicit_sparse_solve_b(typename Active::Tape*,
                                   void* d,
                                   codi::VectorAccessInterface<
                                     typename Active::Real,
                                     typename Active::Identifier>* adj)
{
  using Real = typename Active::Real;

  auto* data = static_cast<ImplicitSparseSolveData<Active>*>(d);
  const int n = data->n;
  const int nnz = data->nnz;

  const size_t maxDim = adj->getVectorSize();

  // Factorize A^T once (passive)
  Eigen::SparseMatrix<double> AT = data->A.transpose();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> luAT;
  luAT.compute(AT);

  for (size_t dim = 0; dim < maxDim; ++dim) {

    // 1) Gather output adjoints (dU_b) and convert to x_b = -(dU_b)
    Eigen::VectorXd x_b(n);
    x_b.setZero();

    for (int i = 0; i < n; ++i) {
      const auto id = data->out_id[i];
      const double dU_b = adj->getAdjoint(id, dim);
      adj->resetAdjoint(id, dim);
      x_b[i] = -dU_b;  // because dU = -x
    }

    // 2) Solve adjoint system: A^T * lambda = x_b
    Eigen::VectorXd lambda = luAT.solve(x_b);

    // 3) RHS adjoint update: b_b += lambda
    for (int i = 0; i < n; ++i) {
      adj->updateAdjoint(data->b_id[i], dim, Real(lambda[i]));
    }

    // 4) Matrix-value adjoint update:
    // Aval_b[k] += -(lambda[row_k] * x[col_k])
    for (int k = 0; k < nnz; ++k) {
      const int r = data->row[k];
      const int c = data->col[k];
      const double contrib = -(lambda[r] * data->x[c]);
      adj->updateAdjoint(data->Aval_id[k], dim, Real(contrib));
    }
  }
};

template<typename Active>
static void implicit_sparse_solve_delete(typename Active::Tape*, void* d)
{
  delete static_cast<ImplicitSparseSolveData<Active>*>(d);
};

// -----------------------------
// Wrapper solve_sys (implicit sparse solve)
// -----------------------------
//template<typename Real>
void solve_sys_sparse(Glob &glob) {
  constexpr int Nsize = 4 * (Ncoords + Nwake);
  using Active = Real;
  using Tape   = typename Active::Tape;

  Tape& tape = Active::getTape();
  const int nnz = glob.R_V_latest;

  // ---- Build passive sparse A from triplets ----
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(nnz);

  for (int k = 0; k < nnz; ++k) {
    triplets.emplace_back(glob.R_V_rows[k],
                          glob.R_V_cols[k],
                          (glob.R_V_vals[k]).getValue());
  }

  Eigen::SparseMatrix<double> A(Nsize, Nsize);
  A.setFromTriplets(triplets.begin(), triplets.end());

  // ---- Build passive rhs b ----
  Eigen::VectorXd b(Nsize);
  for (int i = 0; i < Nsize; ++i) {
    b[i] = (glob.R[i]).getValue();
  }

  // ---- Primal sparse solve: A x = b ----
  Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
  lu.compute(A);

  // You may want error handling here:
  // if (lu.info() != Eigen::Success) { ... }

  Eigen::VectorXd x = lu.solve(b);

  // ---- Write primal output: dU = -x ----
  for (int i = 0; i < Nsize; ++i) {
    glob.dU[i] = Active(-x[i]);
  }

  // ---- If tape inactive: stop here (no extra overhead) ----
  if (!tape.isActive()) return;

  // ---- Tape active: build EF data ----
  auto* data = new ImplicitSparseSolveData<Active>(Nsize, nnz);
  data->A = A;
  data->x = x;

  // Store structure + ids for matrix values
  for (int k = 0; k < nnz; ++k) {
    data->row[k] = glob.R_V_rows[k];
    data->col[k] = glob.R_V_cols[k];
    data->Aval_id[k] = glob.R_V_vals[k].getIdentifier();
  }

  // Store ids for RHS
  for (int i = 0; i < Nsize; ++i) {
    data->b_id[i] = glob.R[i].getIdentifier();
  }

  // Register outputs and store their ids (outputs are glob.dU)
  for (int i = 0; i < Nsize; ++i) {
    tape.registerExternalFunctionOutput(glob.dU[i]);
    data->out_id[i] = glob.dU[i].getIdentifier();
  }

  // Push external function
  tape.pushExternalFunction(
    codi::ExternalFunction<Tape>::create(
      &implicit_sparse_solve_b<Active>,
      data,
      &implicit_sparse_solve_delete<Active>
    )
  );
};
#endif