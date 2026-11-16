// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Tests for the dense matrix helpers and the LDL^T factorization.

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "sna/sna.hpp"
#include "test_util.hpp"

using namespace sna;

namespace {

Matrix<double> randomSymmetric(int n, std::mt19937& rng) {
    std::normal_distribution<double> normal(0.0, 1.0);
    Matrix<double> a(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double value = normal(rng);
            a(i, j) = value;
            a(j, i) = value;
        }
    }
    return a;
}

bool isValidPermutation(const std::vector<int>& perm) {
    std::vector<int> seen(perm.size(), 0);
    for (const int p : perm) {
        if (p < 0 || p >= static_cast<int>(perm.size())) return false;
        if (seen[static_cast<std::size_t>(p)] != 0) return false;
        seen[static_cast<std::size_t>(p)] = 1;
    }
    return true;
}

bool isUnitLower(const Matrix<double>& l, double tolerance) {
    for (int i = 0; i < l.rows(); ++i) {
        for (int j = 0; j < l.cols(); ++j) {
            const double expected = (i == j) ? 1.0 : 0.0;
            if (j > i && std::abs(l(i, j) - expected) > tolerance) return false;
            if (i == j && std::abs(l(i, j) - 1.0) > tolerance) return false;
        }
    }
    return true;
}

bool blocksCoverAll(const LdlResult<double>& factor) {
    int position = 0;
    for (const PivotBlock<double>& block : factor.blocks) {
        if (block.begin != position) return false;
        if (block.size != 1 && block.size != 2) return false;
        position += block.size;
    }
    return position == factor.n;
}

void testMatrixBasics() {
    Matrix<double> a(2, 3);
    a(0, 1) = 4.0;
    CHECK_NEAR(a(0, 1), 4.0, 1e-15);
    CHECK(a.rows() == 2);
    CHECK(a.cols() == 3);

    const Matrix<double> identity = Matrix<double>::identity(4);
    CHECK_NEAR(identity(2, 2), 1.0, 1e-15);
    CHECK_NEAR(identity(2, 3), 0.0, 1e-15);

    const std::vector<double> v = {3.0, 4.0};
    CHECK_NEAR(euclideanNorm(v), 5.0, 1e-15);

    Matrix<double> symmetric(2, 2);
    symmetric(0, 0) = 1.0;
    symmetric(0, 1) = 2.0;
    symmetric(1, 0) = 2.0;
    symmetric(1, 1) = 3.0;
    CHECK(isSymmetric(symmetric));
    symmetric(1, 0) = 2.5;
    CHECK(!isSymmetric(symmetric));
    CHECK(isSymmetric(symmetric, 0.6));
}

void testKnownFactorizations() {
    {
        Matrix<double> a(2, 2);
        a(0, 0) = 2.0;
        a(1, 1) = 2.0;
        a(0, 1) = a(1, 0) = 1.0;
        const LdlResult<double> f = factorSymmetric(a);
        CHECK(reconstructionResidual(a, f) < 1e-14);
        CHECK(isUnitLower(f.L, 1e-15));
        CHECK(blocksCoverAll(f));
        CHECK(isValidPermutation(f.perm));
    }
    {
        // A zero diagonal with a nonzero off-diagonal forces a 2x2 pivot.
        Matrix<double> a(2, 2);
        a(0, 1) = a(1, 0) = 1.0;
        const LdlResult<double> f = factorSymmetric(a);
        CHECK(f.blocks.size() == 1);
        CHECK(f.blocks[0].size == 2);
        CHECK_NEAR(f.blocks[0].a21, 1.0, 1e-15);
        CHECK(reconstructionResidual(a, f) < 1e-15);
    }
    {
        const Matrix<double> zero(3, 3);
        const LdlResult<double> f = factorSymmetric(zero);
        CHECK(f.blocks.size() == 3);
        for (const PivotBlock<double>& block : f.blocks) {
            CHECK(block.size == 1);
            CHECK_NEAR(block.a11, 0.0, 0.0);
        }
        CHECK_NEAR(reconstructionResidual(zero, f), 0.0, 0.0);
    }
    {
        // 3-vertex path at mu = 0. The matrix is singular, and the first column
        // has a zero diagonal with a nonzero entry below it. This case guards
        // the invariant that a zero pivot is chosen only when the remaining
        // block is zero.
        Matrix<double> a(3, 3);
        a(0, 1) = a(1, 0) = 1.0;
        a(0, 2) = a(2, 0) = 1.0;
        const LdlResult<double> f = factorSymmetric(a);
        CHECK(f.blocks.size() == 2);
        CHECK(f.blocks[0].size == 2);
        CHECK_NEAR(f.blocks[1].a11, 0.0, 0.0);
        CHECK(reconstructionResidual(a, f) < 1e-15);
    }
}

void testRandomReconstruction() {
    std::mt19937 rng(12345);
    int two_by_two_blocks = 0;
    double worst = 0.0;

    for (int trial = 0; trial < 300; ++trial) {
        const int n = 1 + static_cast<int>(rng() % 24);
        const Matrix<double> a = randomSymmetric(n, rng);
        const LdlResult<double> f = factorSymmetric(a);

        CHECK(isUnitLower(f.L, 1e-14));
        CHECK(blocksCoverAll(f));
        CHECK(isValidPermutation(f.perm));

        const double residual = reconstructionResidual(a, f);
        worst = std::max(worst, residual);
        CHECK(residual < 1e-11);

        for (const PivotBlock<double>& block : f.blocks) {
            if (block.size == 2) ++two_by_two_blocks;
        }

        // Factorizing the same matrix twice must be deterministic.
        const LdlResult<double> again = factorSymmetric(a);
        CHECK(again.perm == f.perm);
        CHECK(again.blocks.size() == f.blocks.size());
        for (std::size_t i = 0; i < f.blocks.size(); ++i) {
            CHECK(again.blocks[i].begin == f.blocks[i].begin);
            CHECK(again.blocks[i].size == f.blocks[i].size);
            CHECK_NEAR(again.blocks[i].a11, f.blocks[i].a11, 0.0);
            CHECK_NEAR(again.blocks[i].a21, f.blocks[i].a21, 0.0);
            CHECK_NEAR(again.blocks[i].a22, f.blocks[i].a22, 0.0);
        }
    }

    CHECK(two_by_two_blocks > 0);  // the 2x2 branch is covered
    std::cout << "  LDL^T random: worst relative residual " << worst << "\n";
}

void testRankDeficientReconstruction() {
    std::mt19937 rng(999);
    std::normal_distribution<double> normal(0.0, 1.0);
    double worst = 0.0;

    for (int trial = 0; trial < 120; ++trial) {
        const int n = 2 + static_cast<int>(rng() % 18);
        const int rank = 1 + static_cast<int>(rng() % static_cast<unsigned>(n - 1));

        Matrix<double> x(n, rank);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < rank; ++j) x(i, j) = normal(rng);
        }
        Matrix<double> core(rank, rank);
        for (int i = 0; i < rank; ++i) {
            for (int j = 0; j <= i; ++j) {
                const double value = normal(rng);
                core(i, j) = value;
                core(j, i) = value;
            }
        }

        Matrix<double> a(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double sum = 0.0;
                for (int p = 0; p < rank; ++p) {
                    for (int q = 0; q < rank; ++q) sum += x(i, p) * core(p, q) * x(j, q);
                }
                a(i, j) = sum;
            }
        }

        const LdlResult<double> f = factorSymmetric(a);
        const double residual = reconstructionResidual(a, f);
        worst = std::max(worst, residual);
        CHECK(residual < 1e-10);
    }
    std::cout << "  LDL^T rank-deficient: worst relative residual " << worst << "\n";
}

void testTransposedSolve() {
    std::mt19937 rng(4242);
    std::normal_distribution<double> normal(0.0, 1.0);

    for (int trial = 0; trial < 50; ++trial) {
        const int n = 1 + static_cast<int>(rng() % 20);
        Matrix<double> a = randomSymmetric(n, rng);
        for (int i = 0; i < n; ++i) a(i, i) += static_cast<double>(n);  // keep it definite-ish
        const LdlResult<double> f = factorSymmetric(a);

        std::vector<double> rhs(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) rhs[static_cast<std::size_t>(i)] = normal(rng);

        const std::vector<double> x = solveTransposedUnitLower(f.L, rhs);
        double squared_error = 0.0;
        double squared_rhs = 0.0;
        for (int i = 0; i < n; ++i) {
            double lhs = x[static_cast<std::size_t>(i)];
            for (int j = i + 1; j < n; ++j) lhs += f.L(j, i) * x[static_cast<std::size_t>(j)];
            const double diff = lhs - rhs[static_cast<std::size_t>(i)];
            squared_error += diff * diff;
            squared_rhs += rhs[static_cast<std::size_t>(i)] * rhs[static_cast<std::size_t>(i)];
        }
        CHECK(std::sqrt(squared_error) <= 1e-12 * std::max(1.0, std::sqrt(squared_rhs)));
    }
}

void testInvalidInput() {
    CHECK_THROWS(factorSymmetric(Matrix<double>(2, 3)));
}

}  // namespace

void runLdltTests() {
    std::cout << "test_ldlt\n";
    testMatrixBasics();
    testKnownFactorizations();
    testRandomReconstruction();
    testRankDeficientReconstruction();
    testTransposedSolve();
    testInvalidInput();
}
