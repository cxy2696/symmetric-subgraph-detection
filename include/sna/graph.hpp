// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Undirected weighted graph with integer vertex labels. The graph is a
// convenience layer: the detection algorithm operates on a dense symmetric
// matrix, and adjacencyMatrix() produces that matrix in a canonical order
// (labels sorted ascending), so results do not depend on insertion order.

#ifndef SNA_GRAPH_HPP
#define SNA_GRAPH_HPP

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "sna/matrix.hpp"

namespace sna {

template <class T>
class Graph {
public:
    Graph() = default;

    // Create the labels 0, 1, ..., n-1.
    explicit Graph(int n) {
        for (int i = 0; i < n; ++i) addNode(i);
    }

    // Returns true if the label was not present before.
    bool addNode(int label) {
        return adjacency_.emplace(label, std::unordered_map<int, T>()).second;
    }

    // Adding the same edge twice sums the weights. Self-loops are ignored
    // (the label is still created). Missing endpoints are created. An edge with
    // weight 0 is stored and appears in neighbors(), but contributes nothing to
    // adjacencyMatrix() or density(), the two places the detector looks.
    void addEdge(int u, int v, T weight = T(1)) {
        addNode(u);
        if (u == v) return;
        addNode(v);
        adjacency_[u][v] += weight;
        adjacency_[v][u] += weight;
    }

    bool hasNode(int label) const { return adjacency_.count(label) != 0; }

    bool hasEdge(int u, int v) const {
        const auto it = adjacency_.find(u);
        return it != adjacency_.end() && it->second.count(v) != 0;
    }

    T edgeWeight(int u, int v) const {
        const auto it = adjacency_.find(u);
        if (it == adjacency_.end()) return T(0);
        const auto jt = it->second.find(v);
        return jt == it->second.end() ? T(0) : jt->second;
    }

    int numNodes() const { return static_cast<int>(adjacency_.size()); }

    int numEdges() const {
        std::size_t degree_sum = 0;
        for (const auto& entry : adjacency_) degree_sum += entry.second.size();
        return static_cast<int>(degree_sum / 2);
    }

    // Vertex labels in ascending order.
    std::vector<int> labels() const {
        std::vector<int> out;
        out.reserve(adjacency_.size());
        for (const auto& entry : adjacency_) out.push_back(entry.first);
        std::sort(out.begin(), out.end());
        return out;
    }

    std::vector<int> neighbors(int label) const {
        std::vector<int> out;
        const auto it = adjacency_.find(label);
        if (it == adjacency_.end()) return out;
        out.reserve(it->second.size());
        for (const auto& entry : it->second) out.push_back(entry.first);
        std::sort(out.begin(), out.end());
        return out;
    }

    int degree(int label) const {
        const auto it = adjacency_.find(label);
        return it == adjacency_.end() ? 0 : static_cast<int>(it->second.size());
    }

    // Edges (u, v, w) with u < v, sorted by (u, v).
    std::vector<std::tuple<int, int, T>> edges() const {
        std::vector<std::tuple<int, int, T>> out;
        for (const auto& entry : adjacency_) {
            for (const auto& neighbor : entry.second) {
                if (entry.first < neighbor.first) {
                    out.emplace_back(entry.first, neighbor.first, neighbor.second);
                }
            }
        }
        std::sort(out.begin(), out.end());
        return out;
    }

    // Adjacency matrix with rows and columns ordered by labels().
    Matrix<T> adjacencyMatrix() const {
        const std::vector<int> ids = labels();
        const int n = static_cast<int>(ids.size());

        std::unordered_map<int, int> position;
        position.reserve(ids.size() * 2);
        for (int i = 0; i < n; ++i) position[ids[i]] = i;

        Matrix<T> a(n, n);
        for (int i = 0; i < n; ++i) {
            const auto it = adjacency_.find(ids[i]);
            for (const auto& neighbor : it->second) {
                a(i, position.at(neighbor.first)) = neighbor.second;
            }
        }
        return a;
    }

    // Weighted density of the induced subgraph:
    //     2 * sum of internal edge weights / (s * (s - 1)).
    // Equals the usual {0,1} density for unweighted simple graphs. Duplicate
    // labels are removed before the calculation.
    T density(const std::vector<int>& vertices) const {
        std::vector<int> set = vertices;
        std::sort(set.begin(), set.end());
        set.erase(std::unique(set.begin(), set.end()), set.end());

        const int s = static_cast<int>(set.size());
        if (s < 2) return T(0);
        T sum = T(0);
        for (int i = 0; i < s; ++i) {
            const auto it = adjacency_.find(set[static_cast<std::size_t>(i)]);
            if (it == adjacency_.end()) continue;
            for (int j = i + 1; j < s; ++j) {
                const auto jt = it->second.find(set[static_cast<std::size_t>(j)]);
                if (jt != it->second.end()) sum += jt->second;
            }
        }
        return T(2) * sum / (T(s) * T(s - 1));
    }

private:
    std::unordered_map<int, std::unordered_map<int, T>> adjacency_;
};

}  // namespace sna

#endif  // SNA_GRAPH_HPP
