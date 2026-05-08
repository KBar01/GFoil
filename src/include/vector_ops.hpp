#pragma once

// Unified vector/matrix operations for both GFoil_fwd_codi and GFoil_AD.
//
// Template parameter ordering: <int dims..., typename Real>
//   - Allows call sites to specify only the dimension(s) and deduce Real
//     from the pointer arguments: cnp::add_inplace<12>(a, b)
//   - Also accepts explicit Real for call sites that need it:
//     cnp::add_inplace<12, MyReal>(a, b)
//
// equate_block_inplace has no dimension template params so keeps
// template<typename T> only.

#include <cstddef>
#include <type_traits>
#include <cmath>
#ifdef ENABLE_EIGEN
#include <Eigen/Dense>
#endif

namespace cnp {

// ── Scalar ops ───────────────────────────────────────────────────────────────

template<int N, typename Real>
void scalar_mul(const Real* in, Real scalar, Real* out) {
    for (int i = 0; i < N; ++i) out[i] = scalar * in[i];
}

template<int N, typename Real>
void scalar_mul_inplace(Real* in, Real scalar) {
    for (int i = 0; i < N; ++i) in[i] *= scalar;
}

template<int N, typename Real>
void scalar_sub(const Real* in, Real scalar, Real* out) {
    for (int i = 0; i < N; ++i) out[i] = in[i] - scalar;
}

template<int N, typename Real>
void scalar_div(const Real* in, Real scalar, Real* out) {
    for (int i = 0; i < N; ++i) out[i] = in[i] / scalar;
}

template<int N, typename Real>
void scalar_sub_abs(const Real* in, Real scalar, Real* out) {
    for (int i = 0; i < N; ++i) out[i] = std::abs(in[i] - scalar);
}

template<int N, typename Real>
void scalar_sub_abs_inplace(Real* inout, Real scalar) {
    for (int i = 0; i < N; ++i) inout[i] = std::abs(inout[i] - scalar);
}

// ── Concatenation ─────────────────────────────────────────────────────────────

template<int N1, int N2, typename Real>
void hstack(const Real* a, const Real* b, Real* out) {
    for (int i = 0; i < N1; ++i) out[i]      = a[i];
    for (int i = 0; i < N2; ++i) out[N1 + i] = b[i];
}

template<int R1, int R2, int C, typename Real>
void vstack(const Real* A, const Real* B, Real* out) {
    for (int col = 0; col < C; ++col) {
        for (int row = 0; row < R1; ++row)
            out[col*(R1+R2) + row] = A[col*R1 + row];
        for (int row = 0; row < R2; ++row)
            out[col*(R1+R2) + (R1+row)] = B[col*R2 + row];
    }
}

// ── Matrix / vector products ──────────────────────────────────────────────────

template<int M, int N, int P, typename Real>
void matmat_mul(const Real* A, const Real* B, Real* C) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < P; ++j) {
            Real zero = 0.0;
            C[i + j*M] = zero;
            for (int k = 0; k < N; ++k)
                C[i + j*M] += A[i + k*M] * B[k + j*N];
        }
    }
}

template<int M, int N, typename Real>
void outer_product(const Real* a, const Real* b, Real* out) {
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i)
            out[i + j*M] = a[i] * b[j];
}

// ── Vector arithmetic ─────────────────────────────────────────────────────────

template<int N, typename Real>
void add(const Real* a, const Real* b, Real* out) {
    for (int i = 0; i < N; ++i) out[i] = a[i] + b[i];
}

template<int N, typename Real>
void sub(Real* a, const Real* b) {
    for (int i = 0; i < N; ++i) a[i] -= b[i];
}

template<int N, typename Real>
void mul(Real* a, const Real* b, Real* out) {
    for (int i = 0; i < N; ++i) out[i] = a[i] * b[i];
}

template<int N, typename Real>
void add_inplace(Real* a, const Real* b) {
    for (int i = 0; i < N; ++i) a[i] += b[i];
}

template<int N, typename Real>
void equate_inplace(Real* a, const Real* b) {
    for (int i = 0; i < N; ++i) a[i] = b[i];
}

template<int N, typename Real>
void sub_inplace(Real* a, const Real* b) {
    for (int i = 0; i < N; ++i) a[i] -= b[i];
}

template<int N, typename Real>
void mul_inplace(Real* a, const Real* b) {
    for (int i = 0; i < N; ++i) a[i] *= b[i];
}

// ── Block copy ────────────────────────────────────────────────────────────────
// No dimension template params; typename T deduced from arguments.

template<typename T>
void equate_block_inplace(T* A, int ldA, int rowA, int colA,
                          const T* B, int ldB, int rowB, int colB,
                          int nRows, int nCols)
{
    A += colA*ldA + rowA;
    B += colB*ldB + rowB;
    for (int j = 0; j < nCols; ++j) {
        const T* src = B + j*ldB;
        T*       dst = A + j*ldA;
        for (int i = 0; i < nRows; ++i) dst[i] = src[i];
    }
}

}  // namespace cnp
