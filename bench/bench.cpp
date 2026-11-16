// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Scaling check: a disjoint union of 10-cliques, one clique per 10 vertices.
// The dense implementation is O(n^3); this program prints wall-clock time and
// the number of detected groups so changes in cost are visible.

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "sna/sna.hpp"

int main(int argc, char** argv) {
    using namespace sna;
    const int largest = (argc > 1) ? std::atoi(argv[1]) : 800;

    std::cout << "vertices\tseconds\tgroups\n";
    for (int n = 100; n <= largest; n *= 2) {
        Graph<double> graph;
        for (int base = 0; base + 10 <= n; base += 10) {
            for (int i = 0; i < 10; ++i) {
                for (int j = i + 1; j < 10; ++j) graph.addEdge(base + i, base + j);
            }
        }
        // Keep the vertex count exactly n even if n is not a multiple of 10.
        for (int v = 0; v < n; ++v) graph.addNode(v);

        const auto start = std::chrono::steady_clock::now();
        const std::vector<std::vector<int>> groups = detectSymmetries(graph);
        const auto stop = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(stop - start).count();

        std::cout << n << "\t" << seconds << "\t" << groups.size() << "\n";
    }
    return 0;
}
