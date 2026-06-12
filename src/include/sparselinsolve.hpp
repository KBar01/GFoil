#pragma once
#include <cstdint>
#ifdef USE_CODIPACK
#include <codi.hpp>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <vector>
#include <cassert>
#include "real_type.h"
#include "data_structs.h"

// colMajorIndex is inline in real_type.h

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

  // ---- FNV-1a hash of (nnz, rows, cols) for analyzePattern-once cache ----
  constexpr uint64_t FNV_BASIS = 14695981039346656037ULL;
  constexpr uint64_t FNV_PRIME = 1099511628211ULL;
  uint64_t pattern_hash = FNV_BASIS;
  {
    auto mix = [&](int v) {
      const auto* p = reinterpret_cast<const unsigned char*>(&v);
      for (int b = 0; b < 4; ++b) { pattern_hash ^= p[b]; pattern_hash *= FNV_PRIME; }
    };
    mix(nnz);
    for (int k = 0; k < nnz; ++k) { mix(glob.R_V_rows[k]); mix(glob.R_V_cols[k]); }
  }

  // ---- Pattern cache + scatter-add map (Part C) -------------------------
  // All cache state is plain int/bool/uint64_t/double — no Real, no tape
  // interaction. The hash comparison handles both stag-point-move structural
  // changes and restarts from independent solve_coupled calls: a stale cache
  // from a prior solve is safe because a changed pattern triggers re-analyze +
  // slot-map rebuild, and an unchanged pattern means the cached symbolic
  // factorisation and slot map are still valid. Do NOT filter zero values from
  // the triplets — explicit structural zeros keep the pattern constant.
  //
  // On a pattern change we build A with setFromTriplets and record, for each
  // triplet k, the index slot[k] into A.valuePtr() where its (row,col) lives.
  // On an unchanged pattern we skip setFromTriplets and scatter-add (see below)
  // in ascending-k order. Bit-identity vs setFromTriplets was verified over 297
  // solves (golden + slow cases, forward + AD) with a per-iteration memcmp
  // harness; that harness was env-gated and is recoverable from git history.
  static bool     have_pattern = false;
  static uint64_t cached_hash  = 0;
  static int      cached_nnz   = -1;
  static Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
  static Eigen::SparseMatrix<double> A(Nsize, Nsize); // persists across calls
  static std::vector<int>  slot;       // triplet k -> valuePtr index
  static std::vector<char> is_first;   // triplet k is the first to touch its slot

  const bool pattern_changed = !have_pattern
                             || (nnz          != cached_nnz)
                             || (pattern_hash != cached_hash);

  if (pattern_changed) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nnz);
    for (int k = 0; k < nnz; ++k)
      triplets.emplace_back(glob.R_V_rows[k], glob.R_V_cols[k],
                            glob.R_V_vals[k].getValue());
    A.setFromTriplets(triplets.begin(), triplets.end()); // leaves A compressed
    lu.analyzePattern(A);
    cached_hash  = pattern_hash;
    cached_nnz   = nnz;
    have_pattern = true;

    // Build the scatter slot map by walking the compressed CSC structure:
    // inner indices within a column are sorted ascending after makeCompressed,
    // so binary-search each triplet's row within its column's inner range.
    // is_first[k] marks the first (ascending-k) triplet to land on each slot —
    // that triplet ASSIGNS the slot, the rest ADD. This reproduces
    // setFromTriplets exactly: it seeds each cell with the first occurrence's
    // value (preserving the verbatim bits, including -0.0) and then sums any
    // duplicates in triplet order. A plain fill-with-0.0 + accumulate would
    // turn a unique -0.0 entry into +0.0 (0.0 + -0.0 == +0.0), breaking the
    // bitwise match.
    slot.assign(nnz, -1);
    is_first.assign(nnz, 0);
    const int* outer = A.outerIndexPtr();
    const int* inner = A.innerIndexPtr();
    std::vector<char> seen(A.nonZeros(), 0);
    for (int k = 0; k < nnz; ++k) {
      const int c = glob.R_V_cols[k];
      const int r = glob.R_V_rows[k];
      int lo = outer[c], hi = outer[c + 1];
      while (lo < hi) {
        const int mid = lo + ((hi - lo) >> 1);
        const int v = inner[mid];
        if (v == r)     { lo = mid; break; }
        else if (v < r) lo = mid + 1;
        else            hi = mid;
      }
      slot[k] = lo;
      assert(slot[k] < A.nonZeros() && inner[slot[k]] == r);
      if (!seen[slot[k]]) { is_first[k] = 1; seen[slot[k]] = 1; }
    }
  } else {
    // Unchanged pattern: rescatter values into the cached structure. First
    // triplet per slot assigns (seeds verbatim, preserving -0.0); duplicates
    // accumulate in ascending-k order, matching setFromTriplets' summation.
    double* vp = A.valuePtr();
    for (int k = 0; k < nnz; ++k) {
      const double v = glob.R_V_vals[k].getValue();
      if (is_first[k]) vp[slot[k]]  = v;
      else             vp[slot[k]] += v;
    }
  }

  // ---- Build passive rhs b ----
  Eigen::VectorXd b(Nsize);
  for (int i = 0; i < Nsize; ++i) {
    b[i] = (glob.R[i]).getValue();
  }

  // ---- Primal factorize: numerical phase (always called) ----
  lu.factorize(A);

  if (lu.info() != Eigen::Success) {
    for (int i = 0; i < Nsize; ++i)
      glob.dU[i] = Active(0.0);
    return;
  }

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