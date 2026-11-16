// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// End-to-end tests for the shift-and-factorization detector.

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "sna/sna.hpp"
#include "test_util.hpp"

using namespace sna;

namespace {

template <class T>
void addClique(Graph<T>& graph, const std::vector<int>& vertices, T weight = T(1)) {
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        for (std::size_t j = i + 1; j < vertices.size(); ++j) {
            graph.addEdge(vertices[i], vertices[j], weight);
        }
    }
}

Matrix<double> cliqueMatrix(int n, double weight = 1.0) {
    Matrix<double> a(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) a(i, j) = a(j, i) = weight;
    }
    return a;
}

Matrix<double> pathMatrix(int n) {
    Matrix<double> a(n, n);
    for (int i = 0; i + 1 < n; ++i) a(i, i + 1) = a(i + 1, i) = 1.0;
    return a;
}

void testNotebookExample() {
    Graph<double> graph;
    for (int component = 0; component < 5; ++component) {
        std::vector<int> vertices;
        for (int i = 0; i < 10; ++i) vertices.push_back(10 * component + i);
        addClique(graph, vertices);
    }
    // The notebook listed 10 and 19 as extra nodes, but both are already in
    // the cliques; only 50, 90 and 100 are isolated.
    for (const int isolated : {50, 90, 100}) graph.addNode(isolated);

    CHECK(graph.numNodes() == 53);
    CHECK(graph.numEdges() == 225);

    const std::vector<std::vector<int>> expected = {
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
        {10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
        {20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
        {30, 31, 32, 33, 34, 35, 36, 37, 38, 39},
        {40, 41, 42, 43, 44, 45, 46, 47, 48, 49},
    };
    CHECK_GROUPS(detectSymmetries(graph), expected);
}

void testSmallGraphs() {
    {
        Graph<double> graph;
        for (const int base : {0, 3}) {
            addClique(graph, {base, base + 1, base + 2});
        }
        const std::vector<std::vector<int>> expected = {{0, 1, 2}, {3, 4, 5}};
        CHECK_GROUPS(detectSymmetries(graph), expected);
    }
    {
        Graph<double> graph;
        graph.addEdge(7, 9);
        const std::vector<std::vector<int>> expected = {{7, 9}};
        CHECK_GROUPS(detectSymmetries(graph), expected);
    }
    {
        // P4 has no eigenvalue -1 (spectrum 1.618, 0.618, -0.618, -1.618).
        const Matrix<double> path = pathMatrix(4);
        CHECK(detectSymmetries(path).empty());
    }
    {
        // C4 at mu = 0: the null space is spanned by vectors on opposite vertices.
        const Matrix<double> cycle = [&] {
            Matrix<double> a(4, 4);
            for (int i = 0; i < 4; ++i) a(i, (i + 1) % 4) = a((i + 1) % 4, i) = 1.0;
            return a;
        }();
        DetectorParams<double> params;
        params.shift = 0.0;
        const std::vector<std::vector<int>> expected = {{0, 2}, {1, 3}};
        CHECK_GROUPS(detectSymmetries(cycle, params), expected);
    }
}

void testWeightedCliques() {
    {
        Graph<double> graph;
        addClique(graph, {0, 1, 2, 3}, 2.0);
        DetectorParams<double> params;
        params.shift = -2.0;
        const std::vector<std::vector<int>> expected = {{0, 1, 2, 3}};
        CHECK_GROUPS(detectSymmetries(graph, params), expected);
    }
    {
        // With weight 3 the repeated eigenvalue is -3; mu = -1 must not fire.
        Graph<double> graph;
        addClique(graph, {0, 1, 2, 3, 4, 5, 6}, 3.0);
        CHECK(detectSymmetries(graph).empty());
        DetectorParams<double> params;
        params.shift = -3.0;
        const std::vector<std::vector<int>> expected = {{0, 1, 2, 3, 4, 5, 6}};
        CHECK_GROUPS(detectSymmetries(graph, params), expected);
    }
}

void testAttachedClique() {
    // K4 on 0..3; vertex 4 sees all of the clique, vertex 5 only vertex 4.
    // Vectors that sum to zero on 0..3 are null vectors: row 4 sums to zero and
    // row 5 is zero. The clique is symmetric although it is not a component.
    Graph<double> graph;
    addClique(graph, {0, 1, 2, 3});
    for (int i = 0; i < 4; ++i) graph.addEdge(4, i);
    graph.addEdge(5, 4);

    const std::vector<std::vector<int>> expected = {{0, 1, 2, 3}};
    CHECK_GROUPS(detectSymmetries(graph), expected);
}

void testPartialExternalCoupling() {
    // K4 on 0..3. The outside vertex 4 is joined to 0,1 with weight 1 and to
    // 2,3 with weight 2. The full -1 eigenspace of the clique does not survive,
    // but the subspaces on {0,1} and {2,3} do.
    Graph<double> graph;
    addClique(graph, {0, 1, 2, 3});
    graph.addEdge(4, 0, 1.0);
    graph.addEdge(4, 1, 1.0);
    graph.addEdge(4, 2, 2.0);
    graph.addEdge(4, 3, 2.0);

    const std::vector<std::vector<int>> expected = {{0, 1}, {2, 3}};
    CHECK_GROUPS(detectSymmetries(graph), expected);
}

void testNullSpaceBasis() {
    {
        const Matrix<double> clique = cliqueMatrix(5);
        const std::vector<std::vector<double>> basis = approxNullSpace(clique);
        CHECK(basis.size() == 4);  // nullity of K5 + I at mu = -1

        std::vector<int> covered(5, 0);
        for (const std::vector<double>& vector : basis) {
            CHECK_NEAR(euclideanNorm(vector), 1.0, 1e-12);
            CHECK(nullSpaceResidual(clique, -1.0, vector) < 1e-12);
            for (int i = 0; i < 5; ++i) {
                if (std::abs(vector[static_cast<std::size_t>(i)]) > 1e-5) covered[i] = 1;
            }
        }
        for (const int value : covered) CHECK(value == 1);
    }
    {
        // The 2x2 pivot path of P3 at mu = 0 must produce a genuine null vector.
        Matrix<double> path(3, 3);
        path(0, 1) = path(1, 0) = 1.0;
        path(0, 2) = path(2, 0) = 1.0;
        DetectorParams<double> params;
        params.shift = 0.0;

        const std::vector<std::vector<double>> basis = approxNullSpace(path, params);
        CHECK(basis.size() == 1);
        CHECK(nullSpaceResidual(path, 0.0, basis[0]) < 1e-12);

        const std::vector<int> support = [&] {
            std::vector<int> out;
            for (int i = 0; i < 3; ++i) {
                if (std::abs(basis[0][static_cast<std::size_t>(i)]) > 1e-5) out.push_back(i);
            }
            return out;
        }();
        const std::vector<int> expected = {1, 2};
        CHECK(support == expected);
    }
    {
        // A nonsingular shifted matrix has no null vectors.
        const Matrix<double> path = pathMatrix(4);
        CHECK(approxNullSpace(path).empty());
    }
}

void testSupportGrouping() {
    {
        const std::vector<std::vector<double>> vectors = {{1.0, 0.0, 1.0, 0.0},
                                                          {0.0, 1.0, 1.0, 0.0}};
        const std::vector<std::vector<int>> expected = {{0, 1, 2}};
        CHECK_GROUPS(groupBySupport(vectors, 1e-5), expected);
    }
    {
        const std::vector<std::vector<double>> vectors = {{1.0, 0.0, 0.0, 0.0},
                                                          {0.0, 1.0, 0.0, 0.0}};
        CHECK(groupBySupport(vectors, 1e-5).empty());
    }
    {
        const std::vector<std::vector<double>> vectors = {
            {1.0, 0.0, 1e-4, 0.0},
            {0.0, 1.0, 0.0, 1e-4},
        };
        CHECK(groupBySupport(vectors, 1e-3).empty());
        const std::vector<std::vector<int>> expected = {{0, 2}, {1, 3}};
        CHECK_GROUPS(groupBySupport(vectors, 1e-5), expected);
    }
}

void testDensityFilter() {
    Matrix<double> matrix(6, 6);
    for (const int base : {0, 3}) {
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                matrix(base + i, base + j) = matrix(base + j, base + i) = 1.0;
            }
        }
    }

    CHECK_NEAR(inducedDensity(matrix, {0, 1, 2}), 1.0, 1e-15);
    CHECK_NEAR(inducedDensity(matrix, {0, 1, 2, 3}), 0.5, 1e-15);
    CHECK_NEAR(inducedDensity(matrix, {0}), 0.0, 0.0);
    // Duplicate indices must not change the set density.
    CHECK_NEAR(inducedDensity(matrix, {0, 0, 1, 1, 2, 2}), 1.0, 1e-15);

    {
        const std::vector<std::vector<int>> groups = {
            {0, 1, 2}, {0, 1, 2, 3}, {0, 0, 1, 1, 2, 2}};
        const std::vector<std::vector<int>> expected = {{0, 1, 2}};
        CHECK_GROUPS(filterByDensity(matrix, groups, 0.9), expected);
    }
    {
        // A single edge plus an isolated vertex: the set density is 1/3.
        Matrix<double> edge(3, 3);
        edge(0, 1) = edge(1, 0) = 1.0;
        CHECK_NEAR(inducedDensity(edge, {0, 0, 1, 1, 2}), 1.0 / 3.0, 1e-15);
        const std::vector<std::vector<int>> expected = {{0, 1, 2}};
        CHECK_GROUPS(filterByDensity(edge, {{0, 0, 1, 1, 2}}, 0.3), expected);
    }
    {
        Graph<double> graph;
        addClique(graph, {0, 1, 2});
        graph.addNode(9);
        const std::vector<std::vector<int>> groups = {{0, 1, 2, 9}};
        CHECK(filterByDensity(graph, groups, 0.9).empty());  // density 0.5
        CHECK_NEAR(graph.density({0, 1, 2}), 1.0, 1e-15);
        CHECK_NEAR(graph.density({0, 0, 1, 1, 2, 2}), 1.0, 1e-15);
        const std::vector<std::vector<int>> duplicated = {{0, 0, 1, 1, 2, 2}};
        const std::vector<std::vector<int>> expected = {{0, 1, 2}};
        CHECK_GROUPS(filterByDensity(graph, duplicated, 0.9), expected);
    }
}

void testNearZero2x2() {
    // A leading 1x1 pivot sets max|B| = 1, so the relative block tolerance is
    // 1e-8. The trailing 2x2 block is then selected at the second step and has
    // one (first case) or both (second case) eigenvalues below 1e-8.
    DetectorParams<double> params;
    params.shift = 0.0;

    {
        Matrix<double> b(3, 3);
        b(0, 0) = 1.0;
        b(1, 1) = b(2, 2) = 0.7e-8;
        b(1, 2) = b(2, 1) = 1.5e-8;
        const std::vector<std::vector<double>> basis = approxNullSpace(b, params);
        CHECK(basis.size() == 1);
        CHECK(nullSpaceResidual(b, 0.0, basis[0]) < 1e-7);
        const std::vector<std::vector<int>> expected = {{1, 2}};
        CHECK_GROUPS(detectSymmetries(b, params), expected);
    }
    {
        Matrix<double> b(3, 3);
        b(0, 0) = 1.0;
        b(1, 2) = b(2, 1) = 0.9e-8;
        const std::vector<std::vector<double>> basis = approxNullSpace(b, params);
        CHECK(basis.size() == 2);
        for (const std::vector<double>& vector : basis) {
            CHECK(nullSpaceResidual(b, 0.0, vector) < 1e-7);
        }
        const std::vector<std::vector<int>> expected = {{1, 2}};
        CHECK_GROUPS(detectSymmetries(b, params), expected);
    }
}

void testMultiShift() {
    Graph<double> graph;
    addClique(graph, {0, 1, 2});
    for (int i = 0; i < 4; ++i) graph.addEdge(10 + i, 10 + (i + 1) % 4);

    const std::vector<double> shifts = {-1.0, 0.0};
    const std::vector<std::vector<int>> expected = {{0, 1, 2}, {10, 12}, {11, 13}};
    CHECK_GROUPS(detectSymmetries(graph, shifts), expected);

    // Duplicate shifts must not duplicate groups.
    const std::vector<double> repeated = {-1.0, 0.0, -1.0};
    CHECK_GROUPS(detectSymmetries(graph, repeated), expected);
}

void testRandomCliqueUnions() {
    std::mt19937 rng(2026);
    int two_by_two = 0;

    for (int trial = 0; trial < 30; ++trial) {
        const int components = 2 + static_cast<int>(rng() % 4);
        std::vector<std::vector<int>> expected;
        Graph<double> graph;

        int next_label = 0;
        for (int c = 0; c < components; ++c) {
            const int size = 2 + static_cast<int>(rng() % 6);
            std::vector<int> vertices;
            for (int i = 0; i < size; ++i) vertices.push_back(next_label++);
            expected.push_back(vertices);
            addClique(graph, vertices);
        }

        // Relabel by a random permutation so the test also covers label order.
        std::vector<int> permutation(static_cast<std::size_t>(next_label));
        for (int i = 0; i < next_label; ++i) permutation[static_cast<std::size_t>(i)] = i;
        std::shuffle(permutation.begin(), permutation.end(), rng);

        Graph<double> relabeled;
        for (std::size_t c = 0; c < expected.size(); ++c) {
            std::vector<int> vertices;
            for (const int v : expected[c]) {
                vertices.push_back(permutation[static_cast<std::size_t>(v)]);
            }
            addClique(relabeled, vertices);
            std::sort(vertices.begin(), vertices.end());
            expected[c] = vertices;
        }

        const std::vector<std::vector<int>> first = detectSymmetries(relabeled);
        const std::vector<std::vector<int>> second = detectSymmetries(relabeled);
        CHECK_GROUPS(first, expected);
        CHECK_GROUPS(second, expected);
        CHECK(first == second);

        for (const std::vector<int>& group : first) {
            if (group.size() == 2) ++two_by_two;
        }
    }
    CHECK(two_by_two > 0);  // the 2x2-pivot / two-vertex case is covered
}

void testScalarTypes() {
    {
        Graph<float> graph;
        addClique(graph, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
        const std::vector<std::vector<int>> groups = detectSymmetries(graph);
        CHECK(groups.size() == 1);
        CHECK(groups[0].size() == 10);
    }
    {
        Graph<long double> graph;
        addClique(graph, {0, 1, 2, 3, 4, 5});
        const std::vector<std::vector<int>> groups = detectSymmetries(graph);
        CHECK(groups.size() == 1);
        CHECK(groups[0].size() == 6);
    }
}

void testTolerance() {
    // B = J + delta*I has eigenvalues 4 + delta and delta (multiplicity 3).
    // The relative tolerance is block_tol * max|B| = about 4e-8.
    for (const double delta : {1e-10, 1e-6}) {
        Matrix<double> b(4, 4);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) b(i, j) = 1.0;
            b(i, i) += delta;
        }
        DetectorParams<double> params;
        params.shift = 0.0;

        const std::vector<std::vector<double>> basis = approxNullSpace(b, params);
        // max|B| is 4 + delta, so the effective threshold is block_tol * (4 + delta).
        const bool should_detect = delta < params.block_tol * (4.0 + delta);
        if (should_detect) {
            CHECK(basis.size() == 3);
            for (const std::vector<double>& vector : basis) {
                CHECK(nullSpaceResidual(b, 0.0, vector) < 10.0 * delta);
            }
        } else {
            CHECK(basis.empty());
        }
        CHECK(detectSymmetries(b, params).empty() == basis.empty());
    }
}

void testScaleInvariance() {
    // The same structure at three scales must give the same groups: a dense
    // 1e-10 null direction, the same matrix multiplied by 1e6, and a tiny
    // singular matrix.
    DetectorParams<double> params;
    params.shift = 0.0;
    const std::vector<std::vector<int>> expected = {{0, 1, 2, 3}};

    Matrix<double> b(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) b(i, j) = 1.0;
        b(i, i) += 1e-10;
    }
    CHECK_GROUPS(detectSymmetries(b, params), expected);

    Matrix<double> large(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) large(i, j) = 1e6 * b(i, j);
    }
    CHECK_GROUPS(detectSymmetries(large, params), expected);

    Matrix<double> tiny(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) tiny(i, j) = 1e-20;
    }
    CHECK_GROUPS(detectSymmetries(tiny, params), expected);
}

void testDegenerateInputs() {
    CHECK(detectSymmetries(Matrix<double>(0, 0)).empty());
    CHECK(approxNullSpace(Matrix<double>(0, 0)).empty());

    Matrix<double> single(1, 1);
    single(0, 0) = -1.0;
    CHECK(detectSymmetries(single).empty());  // one vertex cannot form a group
    CHECK(approxNullSpace(single).size() == 1);  // but the null vector exists

    // The zero matrix at mu = 0 has a three-dimensional null space. The basis
    // produced by the factorization is the standard basis, so every support is
    // a singleton and no group survives the two-vertex filter. A different
    // basis (e.g. e0 + e1, ...) would give a different partition; this is the
    // basis dependence described in docs/algorithm.md.
    DetectorParams<double> at_zero;
    at_zero.shift = 0.0;
    CHECK(detectSymmetries(Matrix<double>(3, 3), at_zero).empty());

    // A diagonal matrix with one zero entry has a singleton null vector, which
    // does not form a group either.
    Matrix<double> diagonal(3, 3);
    diagonal(0, 0) = 0.0;
    diagonal(1, 1) = 1.0;
    diagonal(2, 2) = 1.0;
    CHECK(detectSymmetries(diagonal, at_zero).empty());

    Graph<double> graph;
    graph.addNode(42);
    CHECK(detectSymmetries(graph).empty());

    CHECK_THROWS(detectSymmetries(Matrix<double>(2, 3)));
    CHECK_THROWS(approxNullSpace(Matrix<double>(2, 3)));
    CHECK_THROWS((Matrix<double>(-1, 2)));

    // Public helpers reject malformed inputs instead of reading out of bounds.
    CHECK_THROWS(groupBySupport(std::vector<std::vector<double>>{{1.0, 0.0}, {0.0}}, 0.5));
    CHECK_THROWS(nullSpaceResidual(Matrix<double>(2, 2), 0.0, std::vector<double>{1.0}));
    CHECK_THROWS(inducedDensity(Matrix<double>(2, 2), std::vector<int>{0, 2}));
    CHECK_THROWS(filterByDensity(Matrix<double>(2, 2), {{0, 2}}, 0.5));
    CHECK_THROWS(reconstructionResidual(Matrix<double>(2, 2), LdlResult<double>{}));
    CHECK_THROWS(solveTransposedUnitLower(Matrix<double>::identity(2), std::vector<double>{1.0}));
    CHECK_THROWS(nullSpaceResidual(Matrix<double>(2, 3), 0.0, std::vector<double>{1.0, 0.0}));
    CHECK_THROWS(inducedDensity(Matrix<double>(2, 3), std::vector<int>{0, 1}));

    // Negative tolerances are rejected instead of silently changing the result.
    DetectorParams<double> bad_params;
    bad_params.block_tol = -1.0;
    CHECK_THROWS(approxNullSpace(Matrix<double>::identity(2), bad_params));
    bad_params.block_tol = 1e-8;
    bad_params.support_tol = -1.0;
    CHECK_THROWS(detectSymmetries(Matrix<double>::identity(2), bad_params));
    CHECK_THROWS(groupBySupport(std::vector<std::vector<double>>{{1.0, 1.0}}, -1.0));
}

}  // namespace

void runDetectorTests() {
    std::cout << "test_detector\n";
    testNotebookExample();
    testSmallGraphs();
    testWeightedCliques();
    testAttachedClique();
    testPartialExternalCoupling();
    testNullSpaceBasis();
    testSupportGrouping();
    testDensityFilter();
    testNearZero2x2();
    testMultiShift();
    testRandomCliqueUnions();
    testScalarTypes();
    testTolerance();
    testScaleInvariance();
    testDegenerateInputs();
}
