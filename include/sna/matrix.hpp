// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Dense row-major matrix. The library only needs a small set of operations:
// construction, element access, and symmetric row/column exchanges used by the
// factorization.

#ifndef SNA_MATRIX_HPP
#define SNA_MATRIX_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace sna {

template <class T>
class Matrix {
public:
    Matrix() = default;

    Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
        if (rows < 0 || cols < 0) {
            throw std::invalid_argument("Matrix: negative dimension");
        }
        data_.assign(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols), T(0));
    }

    static Matrix identity(int n) {
        Matrix m(n, n);
        for (int i = 0; i < n; ++i) m(i, i) = T(1);
        return m;
    }

    int rows() const noexcept { return rows_; }
    int cols() const noexcept { return cols_; }
    bool empty() const noexcept { return rows_ == 0 || cols_ == 0; }

    T& operator()(int i, int j) { return data_[static_cast<std::size_t>(i) * cols_ + j]; }
    const T& operator()(int i, int j) const {
        return data_[static_cast<std::size_t>(i) * cols_ + j];
    }

    T* data() noexcept { return data_.data(); }
    const T* data() const noexcept { return data_.data(); }

    void swapRows(int i, int j) {
        if (i == j || cols_ == 0) return;
        T* a = &(*this)(i, 0);
        T* b = &(*this)(j, 0);
        for (int c = 0; c < cols_; ++c) std::swap(a[c], b[c]);
    }

    void swapCols(int i, int j) {
        if (i == j) return;
        for (int r = 0; r < rows_; ++r) std::swap((*this)(r, i), (*this)(r, j));
    }

private:
    int rows_ = 0;
    int cols_ = 0;
    std::vector<T> data_;
};

// Euclidean norm of a dense vector.
template <class T>
inline T euclideanNorm(const std::vector<T>& v) {
    T sum = T(0);
    for (const T& x : v) sum += x * x;
    return std::sqrt(sum);
}

// Symmetry check with an absolute tolerance. Useful for validating input; the
// factorization and detector assume symmetry and do not check it themselves.
template <class T>
inline bool isSymmetric(const Matrix<T>& a, T tolerance = T(0)) {
    if (a.rows() != a.cols()) return false;
    for (int i = 0; i < a.rows(); ++i) {
        for (int j = i + 1; j < a.cols(); ++j) {
            if (std::abs(a(i, j) - a(j, i)) > tolerance) return false;
        }
    }
    return true;
}

}  // namespace sna

#endif  // SNA_MATRIX_HPP
