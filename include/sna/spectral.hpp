// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Detection of vertex groups that support vectors of an eigenspace of a graph at
// a chosen shift (Algorithm 1 of the accompanying note):
//
//   lines 4-5    B = A - mu*I
//   line 6       P B P^T = L*D*L^T
//   lines 7-35   near-zero blocks of D -> seeds -> approximate null vectors
//   lines 37-43  union of supports -> connected components
//   lines 44-50  refinement (here: a density filter, applied by the caller)
//
// The pipeline is exposed in three pieces:
//   approxNullSpace()  factorization + seed recovery
//   groupBySupport()   supports + connected components
//   detectSymmetries() the two composed
//
// Matrix overloads use 0-based matrix indices. Graph overloads translate the
// result back to the graph's vertex labels.

#ifndef SNA_SPECTRAL_HPP
#define SNA_SPECTRAL_HPP

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "sna/graph.hpp"
#include "sna/ldlt.hpp"
#include "sna/matrix.hpp"

namespace sna {

template <class T>
struct DetectorParams {
    T shift = T(-1);          // mu: candidate eigenvalue
    T block_tol = T(1e-8);    // epsilon, relative to max|B|: a smaller pivot/eigenvalue is zero
    T support_tol = T(1e-5);  // tau: support threshold on a normalized null vector
};

namespace detail {

class DisjointSets {
public:
    explicit DisjointSets(int n)
        : parent_(static_cast<std::size_t>(n)), size_(static_cast<std::size_t>(n), 1) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    int find(int x) {
        while (parent_[static_cast<std::size_t>(x)] != x) {
            parent_[static_cast<std::size_t>(x)] =
                parent_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(x)])];
            x = parent_[static_cast<std::size_t>(x)];
        }
        return x;
    }

    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (size_[static_cast<std::size_t>(a)] < size_[static_cast<std::size_t>(b)]) {
            std::swap(a, b);
        }
        parent_[static_cast<std::size_t>(b)] = a;
        size_[static_cast<std::size_t>(a)] += size_[static_cast<std::size_t>(b)];
    }

private:
    std::vector<int> parent_;
    std::vector<int> size_;
};

// Sorted, duplicate-free copy of a vertex list.
inline std::vector<int> normalizedGroup(const std::vector<int>& group) {
    std::vector<int> out = group;
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

// Density of a duplicate-free vertex set; the caller has validated the indices.
template <class T>
T densityOfSet(const Matrix<T>& a, const std::vector<int>& vertices) {
    const int s = static_cast<int>(vertices.size());
    if (s < 2) return T(0);
    T sum = T(0);
    for (int i = 0; i < s; ++i) {
        for (int j = i + 1; j < s; ++j) {
            sum += a(vertices[static_cast<std::size_t>(i)],
                     vertices[static_cast<std::size_t>(j)]);
        }
    }
    return T(2) * sum / (T(s) * T(s - 1));
}

}  // namespace detail

// A basis of the approximate null space of B = A - shift*I, in original vertex
// order. Every vector has unit Euclidean norm. One vector is produced per
// near-zero 1x1 block and per near-zero eigenvalue of a 2x2 block.
template <class T>
std::vector<std::vector<T>> approxNullSpace(const Matrix<T>& a,
                                            const DetectorParams<T>& params = {}) {
    static_assert(std::is_floating_point<T>::value,
                  "sna requires a floating-point scalar type");
    if (a.rows() != a.cols()) {
        throw std::invalid_argument("approxNullSpace: matrix must be square");
    }
    if (!(params.block_tol >= T(0)) || !(params.support_tol >= T(0))) {
        throw std::invalid_argument("approxNullSpace: tolerances must be non-negative");
    }

    const int n = a.rows();
    std::vector<std::vector<T>> basis;
    if (n == 0) return basis;

    // Algorithm 1, lines 4-6.
    Matrix<T> shifted = a;
    for (int i = 0; i < n; ++i) shifted(i, i) -= params.shift;
    const LdlResult<T> factor = factorSymmetric(shifted);

    // A pivot is zero relative to the scale of the matrix that was factored.
    // The paper uses an absolute epsilon; making it relative keeps the test
    // meaningful for weighted or rescaled matrices (max|B| = 1 for 0/1 inputs
    // with a shift of order one, which is the case the paper discusses).
    const T matrix_scale = maxAbsEntry(shifted);
    const T block_tol = params.block_tol * (matrix_scale > T(0) ? matrix_scale : T(1));

    std::vector<T> seed(static_cast<std::size_t>(n), T(0));

    // Algorithm 1, lines 30-35: solve L^T v = u, normalize, undo the permutation.
    auto addSeed = [&](int position, T x, T y) {
        std::fill(seed.begin(), seed.end(), T(0));
        seed[static_cast<std::size_t>(position)] = x;
        if (y != T(0)) seed[static_cast<std::size_t>(position + 1)] = y;

        std::vector<T> v = solveTransposedUnitLower(factor.L, seed);
        const T scale = euclideanNorm(v);
        if (!(scale > T(0))) return;

        std::vector<T> vector(static_cast<std::size_t>(n), T(0));
        for (int i = 0; i < n; ++i) {
            vector[static_cast<std::size_t>(factor.perm[static_cast<std::size_t>(i)])] =
                v[static_cast<std::size_t>(i)] / scale;
        }
        basis.push_back(std::move(vector));
    };

    // Algorithm 1, lines 7-28. A 1x1 block contributes e_k. A 2x2 block
    // contributes its eigenvector when an eigenvalue is near zero; using e_k
    // would not give a null vector in general.
    for (const PivotBlock<T>& block : factor.blocks) {
        if (block.size == 1) {
            if (std::abs(block.a11) < block_tol) addSeed(block.begin, T(1), T(0));
            continue;
        }

        const T a11 = block.a11;
        const T a21 = block.a21;
        const T a22 = block.a22;
        const T trace = a11 + a22;
        const T disc = std::hypot(a11 - a22, T(2) * a21);
        const T eigenvalues[2] = {(trace + disc) / T(2), (trace - disc) / T(2)};

        for (const T eigenvalue : eigenvalues) {
            if (std::abs(eigenvalue) >= block_tol) continue;

            // Eigenvector of [a11 a21; a21 a22]. The two expressions are
            // equivalent; pick the one that does not subtract close values.
            T x = T(0);
            T y = T(0);
            if (std::abs(a21) >= std::abs(eigenvalue - a11)) {
                x = a21;
                y = eigenvalue - a11;
            } else {
                x = eigenvalue - a22;
                y = a21;
            }
            const T length = std::hypot(x, y);
            if (length == T(0)) continue;  // impossible for a 2x2 pivot with a21 != 0
            addSeed(block.begin, x / length, y / length);
        }
    }
    return basis;
}

// Algorithm 1, lines 37-43. Each vector contributes a clique over its support;
// the connected components of that support graph are returned. A star over the
// support spans the same connectivity, so only |supp| unions are needed.
//
// All vectors must have the same length (the vertex count). The vectors are used
// as given; approxNullSpace returns unit vectors, which is what makes
// support_tol a relative threshold.
template <class T>
std::vector<std::vector<int>> groupBySupport(const std::vector<std::vector<T>>& vectors,
                                             T support_tol) {
    if (!(support_tol >= T(0))) {
        throw std::invalid_argument("groupBySupport: support_tol must be non-negative");
    }
    if (vectors.empty()) return {};

    const std::size_t length = vectors.front().size();
    for (const std::vector<T>& vector : vectors) {
        if (vector.size() != length) {
            throw std::invalid_argument("groupBySupport: vectors must have equal length");
        }
    }

    const int n = static_cast<int>(length);
    detail::DisjointSets sets(n);

    for (const std::vector<T>& vector : vectors) {
        int first = -1;
        for (int i = 0; i < n; ++i) {
            if (std::abs(vector[static_cast<std::size_t>(i)]) > support_tol) {
                if (first < 0) {
                    first = i;
                } else {
                    sets.unite(first, i);
                }
            }
        }
    }

    std::vector<std::vector<int>> buckets(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) buckets[static_cast<std::size_t>(sets.find(i))].push_back(i);

    std::vector<std::vector<int>> groups;
    for (std::vector<int>& bucket : buckets) {
        if (bucket.size() >= 2) groups.push_back(std::move(bucket));
    }
    std::sort(groups.begin(), groups.end());
    return groups;
}

// Algorithm 1 on a dense symmetric matrix. Returns disjoint groups of at least
// two indices, each sorted ascending; the list is sorted lexicographically.
// The matrix is assumed symmetric and is not checked for symmetry.
template <class T>
std::vector<std::vector<int>> detectSymmetries(const Matrix<T>& a,
                                              const DetectorParams<T>& params = {}) {
    return groupBySupport(approxNullSpace(a, params), params.support_tol);
}

// Run the detector for several shifts and merge the results. Groups found at
// more than one shift are returned once (Algorithm 1, line 51); the result is a
// union of per-shift groups, not necessarily a partition.
template <class T>
std::vector<std::vector<int>> detectSymmetries(const Matrix<T>& a,
                                              const std::vector<T>& shifts,
                                              DetectorParams<T> params = {}) {
    std::set<std::vector<int>> unique;
    for (const T shift : shifts) {
        params.shift = shift;
        const std::vector<std::vector<int>> found = detectSymmetries(a, params);
        unique.insert(found.begin(), found.end());
    }
    return std::vector<std::vector<int>>(unique.begin(), unique.end());
}

// ||(A - shift*I) x|| / (||x|| * max|A - shift*I|). The denominator uses the
// same matrix scale as block_tol, so a "small" value means the same thing at
// any scale. Returns 0 for x = 0, where the quotient is undefined.
template <class T, class U>
T nullSpaceResidual(const Matrix<T>& a, U shift, const std::vector<T>& x) {
    if (a.rows() != a.cols()) {
        throw std::invalid_argument("nullSpaceResidual: matrix must be square");
    }
    const int n = a.rows();
    if (static_cast<int>(x.size()) != n) {
        throw std::invalid_argument("nullSpaceResidual: vector size does not match matrix");
    }
    if (n == 0) return T(0);
    const T vector_scale = euclideanNorm(x);
    if (!(vector_scale > T(0))) return T(0);
    const T mu = static_cast<T>(shift);

    T sum = T(0);
    T matrix_scale = T(0);
    for (int i = 0; i < n; ++i) {
        T row = -mu * x[static_cast<std::size_t>(i)];
        for (int j = 0; j < n; ++j) {
            row += a(i, j) * x[static_cast<std::size_t>(j)];
            const T entry = a(i, j) - (i == j ? mu : T(0));
            matrix_scale = std::max(matrix_scale, std::abs(entry));
        }
        sum += row * row;
    }
    const T denominator = vector_scale * (matrix_scale > T(0) ? matrix_scale : T(1));
    return std::sqrt(sum) / denominator;
}

// Algorithm 1 on a graph. Returns the original vertex labels.
template <class T>
std::vector<std::vector<int>> detectSymmetries(const Graph<T>& graph,
                                              const DetectorParams<T>& params = {}) {
    const std::vector<int> labels = graph.labels();
    const std::vector<std::vector<int>> index_groups =
        detectSymmetries(graph.adjacencyMatrix(), params);

    std::vector<std::vector<int>> out;
    out.reserve(index_groups.size());
    for (const std::vector<int>& group : index_groups) {
        std::vector<int> vertices;
        vertices.reserve(group.size());
        for (const int index : group) vertices.push_back(labels[static_cast<std::size_t>(index)]);
        out.push_back(std::move(vertices));
    }
    return out;
}

template <class T>
std::vector<std::vector<int>> detectSymmetries(const Graph<T>& graph,
                                              const std::vector<T>& shifts,
                                              DetectorParams<T> params = {}) {
    const std::vector<int> labels = graph.labels();
    const Matrix<T> adjacency = graph.adjacencyMatrix();

    std::set<std::vector<int>> unique;
    for (const T shift : shifts) {
        params.shift = shift;
        for (const std::vector<int>& group : detectSymmetries(adjacency, params)) {
            std::vector<int> vertices;
            vertices.reserve(group.size());
            for (const int index : group) {
                vertices.push_back(labels[static_cast<std::size_t>(index)]);
            }
            unique.insert(std::move(vertices));
        }
    }
    return std::vector<std::vector<int>>(unique.begin(), unique.end());
}

// Weighted density of the submatrix induced by `vertices`:
//     2 * sum of off-diagonal entries / (s*(s-1)).
// Duplicate indices are removed before the calculation; indices must be in
// [0, a.rows()).
template <class T>
T inducedDensity(const Matrix<T>& a, const std::vector<int>& vertices) {
    if (a.rows() != a.cols()) {
        throw std::invalid_argument("inducedDensity: matrix must be square");
    }
    const std::vector<int> set = detail::normalizedGroup(vertices);
    for (const int vertex : set) {
        if (vertex < 0 || vertex >= a.rows()) {
            throw std::invalid_argument("inducedDensity: vertex index out of range");
        }
    }
    return detail::densityOfSet(a, set);
}

// Keep the groups (after removing duplicate vertices) whose weighted induced
// density is at least min_density. This filters; it does not split a group.
template <class T>
std::vector<std::vector<int>> filterByDensity(const Matrix<T>& a,
                                              const std::vector<std::vector<int>>& groups,
                                              T min_density) {
    if (a.rows() != a.cols()) {
        throw std::invalid_argument("filterByDensity: matrix must be square");
    }
    std::vector<std::vector<int>> out;
    for (const std::vector<int>& group : groups) {
        const std::vector<int> set = detail::normalizedGroup(group);
        if (set.size() < 2) continue;
        for (const int vertex : set) {
            if (vertex < 0 || vertex >= a.rows()) {
                throw std::invalid_argument("filterByDensity: vertex index out of range");
            }
        }
        if (detail::densityOfSet(a, set) < min_density) continue;
        out.push_back(set);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

// Graph variant of the density filter; groups contain vertex labels.
template <class T>
std::vector<std::vector<int>> filterByDensity(const Graph<T>& graph,
                                              const std::vector<std::vector<int>>& groups,
                                              T min_density) {
    std::vector<std::vector<int>> out;
    for (const std::vector<int>& group : groups) {
        const std::vector<int> set = detail::normalizedGroup(group);
        if (set.size() < 2) continue;
        if (graph.density(set) < min_density) continue;
        out.push_back(set);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

}  // namespace sna

#endif  // SNA_SPECTRAL_HPP
