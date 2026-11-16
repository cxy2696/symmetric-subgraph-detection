// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// LDL^T factorization of a real symmetric matrix, with 1x1 and 2x2 pivots:
//
//     A(perm[i], perm[j]) = (L * D * L^T)(i, j).
//
// L is unit lower triangular. D is block diagonal with 1x1 and 2x2 blocks and is
// stored as a list of blocks, which is all the symmetry detector needs.
//
// Pivot strategy
// --------------
// B = A - mu*I is singular by construction, so the factorization has to be able
// to emit a zero pivot. That is safe only when the remaining row and column are
// zero; otherwise elimination would divide by zero. Complete pivoting over the
// trailing block (Bunch--Parlett), modified to prefer diagonal pivots, gives
// this invariant directly:
//
//     let M = max |B[i][j]| (i != j) and G = max |B[i][i]| over the trailing block.
//     If G == 0 and M == 0: the trailing block is zero; emit zero 1x1 pivots.
//     If G >= 0.5*M:         use a 1x1 pivot on the largest diagonal entry.
//     Otherwise:             use the 2x2 block holding the largest off-diagonal entry.
//
// Consequences, all used by the proof in docs/algorithm.md:
//   * Every 1x1 pivot used for elimination is nonzero (it is G or at least 0.5*M).
//   * Every 2x2 pivot has |det| >= (1 - 0.25)*M^2 = 0.75*M^2, so it is invertible
//     and the multipliers are bounded by 2.
//   * In exact arithmetic, a zero pivot appears only when every entry of the
//     trailing block is zero, so no step divides by zero.
//
// Complete pivoting costs O(n^3) comparisons in addition to the O(n^3)
// elimination; the asymptotic complexity is unchanged.

#ifndef SNA_LDLT_HPP
#define SNA_LDLT_HPP

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "sna/matrix.hpp"

namespace sna {

// One diagonal block of D. For a 1x1 block only a11 is used. For a 2x2 block
// the entries are [a11 a21; a21 a22].
template <class T>
struct PivotBlock {
    int begin;  // position of the block in the permuted ordering
    int size;   // 1 or 2
    T a11;
    T a21;
    T a22;
};

template <class T>
struct LdlResult {
    Matrix<T> L;                    // unit lower triangular
    std::vector<PivotBlock<T>> blocks;
    std::vector<int> perm;          // perm[i] = original index at position i
    int n = 0;
};

// Factor a symmetric matrix. Throws std::invalid_argument if A is not square.
template <class T>
LdlResult<T> factorSymmetric(const Matrix<T>& A) {
    static_assert(std::is_floating_point<T>::value,
                  "sna requires a floating-point scalar type");

    if (A.rows() != A.cols()) {
        throw std::invalid_argument("factorSymmetric: matrix must be square");
    }

    const int n = A.rows();
    Matrix<T> S = A;                 // trailing Schur complement, kept symmetric
    Matrix<T> L = Matrix<T>::identity(n);
    std::vector<int> perm(static_cast<std::size_t>(n));
    std::iota(perm.begin(), perm.end(), 0);
    std::vector<PivotBlock<T>> blocks;
    blocks.reserve(static_cast<std::size_t>(n));

    // Symmetric interchange of positions i and j. The first `done` columns of L
    // hold multipliers that are already final; only the corresponding entries of
    // rows i and j are exchanged, so L stays lower triangular with unit diagonal.
    auto swapPositions = [&](int i, int j, int done) {
        if (i == j) return;
        S.swapRows(i, j);
        S.swapCols(i, j);
        for (int c = 0; c < done; ++c) std::swap(L(i, c), L(j, c));
        std::swap(perm[static_cast<std::size_t>(i)], perm[static_cast<std::size_t>(j)]);
    };

    const T diagonal_preference = T(0.5);

    int k = 0;
    while (k < n) {
        T max_diag = T(0);
        int diag_pos = -1;
        T max_off = T(0);
        int off_i = -1;
        int off_j = -1;

        for (int i = k; i < n; ++i) {
            const T ai = std::abs(S(i, i));
            if (ai > max_diag) {
                max_diag = ai;
                diag_pos = i;
            }
            for (int j = i + 1; j < n; ++j) {
                const T aj = std::abs(S(i, j));
                if (aj > max_off) {
                    max_off = aj;
                    off_i = i;
                    off_j = j;
                }
            }
        }

        if (max_diag == T(0) && max_off == T(0)) {
            // The trailing block is exactly zero. Every remaining position is a
            // (legitimate) zero 1x1 pivot; no elimination is needed.
            for (; k < n; ++k) blocks.push_back(PivotBlock<T>{k, 1, T(0), T(0), T(0)});
            break;
        }

        const bool use_1x1 =
            max_diag > T(0) && (max_off == T(0) || max_diag >= diagonal_preference * max_off);

        if (use_1x1) {
            swapPositions(k, diag_pos, k);
            const T d = S(k, k);  // nonzero by the pivot rule
            blocks.push_back(PivotBlock<T>{k, 1, d, T(0), T(0)});

            for (int i = k + 1; i < n; ++i) L(i, k) = S(i, k) / d;
            for (int i = k + 1; i < n; ++i) {
                for (int j = k + 1; j <= i; ++j) {
                    S(i, j) -= L(i, k) * S(j, k);
                    S(j, i) = S(i, j);
                }
            }
            ++k;
        } else {
            // Move the pair (off_i, off_j) to positions (k, k+1).
            swapPositions(k, off_i, k);
            swapPositions(k + 1, off_j, k);

            const T a = S(k, k);
            const T b = S(k + 1, k);
            const T c = S(k + 1, k + 1);
            const T det = a * c - b * b;  // |det| >= 0.75 * b^2 > 0
            blocks.push_back(PivotBlock<T>{k, 2, a, b, c});

            for (int i = k + 2; i < n; ++i) {
                const T x = S(i, k);
                const T y = S(i, k + 1);
                L(i, k) = (c * x - b * y) / det;
                L(i, k + 1) = (-b * x + a * y) / det;
            }
            for (int i = k + 2; i < n; ++i) {
                const T l1 = L(i, k);
                const T l2 = L(i, k + 1);
                for (int j = k + 2; j <= i; ++j) {
                    const T m1 = L(j, k);
                    const T m2 = L(j, k + 1);
                    S(i, j) -= l1 * (a * m1 + b * m2) + l2 * (b * m1 + c * m2);
                    S(j, i) = S(i, j);
                }
            }
            k += 2;
        }
    }

    return LdlResult<T>{std::move(L), std::move(blocks), std::move(perm), n};
}

// Solve L^T x = rhs, where L is unit lower triangular. Back substitution from
// the last row; the unit diagonal needs no division.
template <class T>
std::vector<T> solveTransposedUnitLower(const Matrix<T>& L, const std::vector<T>& rhs) {
    if (L.rows() != L.cols()) {
        throw std::invalid_argument("solveTransposedUnitLower: matrix must be square");
    }
    if (static_cast<int>(rhs.size()) != L.rows()) {
        throw std::invalid_argument("solveTransposedUnitLower: right-hand side size mismatch");
    }
    const int n = L.rows();
    std::vector<T> x(static_cast<std::size_t>(n), T(0));
    for (int i = n - 1; i >= 0; --i) {
        T acc = rhs[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < n; ++j) acc -= L(j, i) * x[static_cast<std::size_t>(j)];
        x[static_cast<std::size_t>(i)] = acc;
    }
    return x;
}

// Relative Frobenius residual
//     || A(perm[i],perm[j]) - (L*D*L^T)(i,j) ||_F / max(1, ||A||_F).
// O(n^3) time and O(n^2) workspace; intended for validation and debugging.
template <class T>
T reconstructionResidual(const Matrix<T>& A, const LdlResult<T>& f) {
    if (A.rows() != A.cols() || f.n != A.rows() || f.L.rows() != f.n || f.L.cols() != f.n ||
        static_cast<int>(f.perm.size()) != f.n) {
        throw std::invalid_argument("reconstructionResidual: incompatible matrix and factor");
    }
    int covered = 0;
    for (const PivotBlock<T>& block : f.blocks) {
        if (block.begin != covered || (block.size != 1 && block.size != 2)) {
            throw std::invalid_argument("reconstructionResidual: malformed pivot blocks");
        }
        covered += block.size;
    }
    if (covered != f.n) {
        throw std::invalid_argument("reconstructionResidual: pivot blocks do not cover the matrix");
    }

    const int n = A.rows();
    if (n == 0) return T(0);

    // T1 = D * L^T. D is block diagonal, so this costs O(n^2).
    Matrix<T> t1(n, n);
    for (const PivotBlock<T>& blk : f.blocks) {
        const int p = blk.begin;
        if (blk.size == 1) {
            for (int j = 0; j < n; ++j) t1(p, j) = blk.a11 * f.L(j, p);
        } else {
            for (int j = 0; j < n; ++j) {
                t1(p, j) = blk.a11 * f.L(j, p) + blk.a21 * f.L(j, p + 1);
                t1(p + 1, j) = blk.a21 * f.L(j, p) + blk.a22 * f.L(j, p + 1);
            }
        }
    }

    T num = T(0);
    T den = T(0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            T acc = T(0);
            for (int p = 0; p <= i; ++p) acc += f.L(i, p) * t1(p, j);
            const T ref = A(f.perm[static_cast<std::size_t>(i)],
                            f.perm[static_cast<std::size_t>(j)]);
            const T diff = acc - ref;
            num += diff * diff;
            den += ref * ref;
        }
    }
    return std::sqrt(num) / std::max(T(1), std::sqrt(den));
}

}  // namespace sna

#endif  // SNA_LDLT_HPP
