// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Reference example: five disjoint 10-cliques plus isolated vertices, shifted by
// mu = -1. This is the graph from the original notebook; the detector should
// return the five cliques and nothing else.

#include <iostream>

#include "sna/sna.hpp"

int main() {
    using namespace sna;

    Graph<double> graph;
    for (int component = 0; component < 5; ++component) {
        const int base = 10 * component;
        for (int i = 0; i < 10; ++i) {
            for (int j = i + 1; j < 10; ++j) graph.addEdge(base + i, base + j);
        }
    }
    // The notebook also listed 10 and 19 as extra nodes, but both are already
    // in the cliques; only 50, 90 and 100 are isolated.
    for (const int isolated : {50, 90, 100}) graph.addNode(isolated);

    std::cout << "vertices: " << graph.numNodes() << ", edges: " << graph.numEdges() << "\n";

    const std::vector<std::vector<int>> groups = detectSymmetries(graph);
    std::cout << "groups: " << groups.size() << "\n";
    for (const std::vector<int>& group : groups) {
        std::cout << "  " << group.size() << " vertices:";
        for (const int vertex : group) std::cout << " " << vertex;
        std::cout << "\n";
    }

    // Audit: factor the shifted adjacency matrix and report the reconstruction
    // residual. This checks the numerical core, not the detection result.
    Matrix<double> shifted = graph.adjacencyMatrix();
    const DetectorParams<double> params;
    for (int i = 0; i < shifted.rows(); ++i) shifted(i, i) -= params.shift;
    const LdlResult<double> factor = factorSymmetric(shifted);
    std::cout << "LDL^T relative residual: " << reconstructionResidual(shifted, factor) << "\n";
    return 0;
}
