#pragma once

#include <cstddef>  // for size_t
#include <type_traits>
#include "real_type.hpp"
#ifdef ENABLE_EIGEN
#include <Eigen/Dense>
#endif

namespace cnp {


// ========================
//  Scalar Ops
// ========================

template<typename Real, int N>
void scalar_mul(const Real* in, Real scalar, Real* out) {
    for (int i = 0; i < N; ++i)
        out[i] = scalar * in[i];
}


template<typename Real, int N>
void scalar_mul_inplace(Real* in, Real scalar) {
    for (int i = 0; i < N; ++i)
        in[i] *= scalar;
}

template<typename Real, int N>
void scalar_sub(const Real* in, Real scalar, Real* out) {
    for (int i = 0; i < N; ++i)
        out[i] = in[i] - scalar;
}

template<typename Real, int N>
void scalar_div(const Real* in, Real scalar, Real* out) {
    for (int i = 0; i < N; ++i)
        out[i] = in[i] / scalar;
}


template<typename Real, int N>
void scalar_sub_abs(const Real* in, Real scalar, Real* out) {
    for (int i = 0; i < N; ++i)
        out[i] = std::abs(in[i] - scalar);
}

template<typename Real, int N>
void scalar_sub_abs_inplace(Real* inout, Real scalar) {
    for (int i = 0; i < N; ++i)
        inout[i] = std::abs(inout[i] - scalar);
}

// ========================
//  Concatenation
// ========================

template<typename Real, int N1, int N2>
void hstack(const Real* a, const Real* b, Real* out) {
    for (int i = 0; i < N1; ++i)
        out[i] = a[i];
    for (int i = 0; i < N2; ++i)
        out[N1 + i] = b[i];
}


template<typename Real, int R1, int R2, int C>
void vstack(const Real* A, const Real* B, Real* out) {
    for (int col = 0; col < C; ++col) {
        for (int row = 0; row < R1; ++row)
            out[col * (R1 + R2) + row] = A[col * R1 + row];

        for (int row = 0; row < R2; ++row)
            out[col * (R1 + R2) + (R1 + row)] = B[col * R2 + row];
    }
}

// =========================
// Manual dot product
// =========================

template<typename Real, int M, int N, int P>
void matmat_mul(const Real* A, const Real* B, Real* C) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < P; ++j) {
            Real zero = 0.0;
            C[i + j*M] = zero;
            for (int k = 0; k < N; ++k) {
                C[i + j*M] += A[i + k*M] * B[k + j*N];
            }
        }
    }
}

// Outer product: C = a * b^T, size(C) = M x N
template<typename Real, int M, int N>
void outer_product(const Real* a, const Real* b, Real* out) {
    // Outer product: C = a * b^T, size(C) = M x N
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            out[i + j*M] = a[i] * b[j];
        }
    }
}


// Vector addition: out = a + b
template<typename Real, int N>
void add(const Real* a, const Real* b, Real* out) {
    for (int i = 0; i < N; ++i) {
        out[i] = a[i] + b[i];
    }
}

template<typename Real, int N>
void sub(Real* a, const Real* b) {
    for (int i = 0; i < N; ++i)
        a[i] -= b[i];
}

template<typename Real, int N>
void mul(Real* a, const Real* b, Real* out) {
    for (int i = 0; i < N; ++i)
        out[i] = a[i] * b[i];
}

template<typename Real, int N>
void add_inplace(Real* a, const Real* b) {
    // In-place vector addition: a += b
    for (int i = 0; i < N; ++i) {
        a[i] += b[i];
    }
}

template<typename Real, int N>
void equate_inplace(Real* a, const Real* b) {
    // In-place vector assigning: a = b
    for (int i = 0; i < N; ++i) {
        a[i] = b[i];
    }
}

/// Add a rectangular block (nRows × nCols, column‑major) from B into A.
///
/// * **A, B** – pointers to the *first* element of each matrix  
/// * **ldA, ldB** – leading dimensions (the physical #rows – the stride from
///                  one element of column k to column k+1)  
/// * **rowA, colA** – top‑left corner of the destination block inside A  
/// * **rowB, colB** – top‑left corner of the source     block inside B  
/// * **nRows, nCols** – size of the block to add
///
/// All arrays are assumed **column‑major**.

template <typename Real>
void equate_block_inplace(Real* A, int ldA, int rowA, int colA,
                    const Real* B, int ldB, int rowB, int colB,
                    int nRows, int nCols){

    // Move pointers to the start of each block
    A += colA * ldA + rowA;
    B += colB * ldB + rowB;

    // Walk column by column (the contiguous direction)
    for (int j = 0; j < nCols; ++j){
        
        const Real* src = B + j * ldB;
        Real*       dst = A + j * ldA;

        for (int i = 0; i < nRows; ++i){dst[i] = src[i];}
    }
}


template<typename Real, int N>
void sub_inplace(Real* a, const Real* b) {
    for (int i = 0; i < N; ++i)
        a[i] -= b[i];
}

template<typename Real, int N>
void mul_inplace(Real* a, const Real* b) {
    for (int i = 0; i < N; ++i)
        a[i] *= b[i];
}


}  // namespace vecops
