// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Command-line driver for an edge list:
//
//     edge_list <file> [shift]
//
// Each line is "u v" or "u v w"; '#' starts a comment. Vertices that appear in
// no edge cannot be represented in this format and are not added to the graph.

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sna/sna.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <edge-list> [shift]\n";
        return 2;
    }

    std::ifstream input(argv[1]);
    if (!input) {
        std::cerr << "cannot open " << argv[1] << "\n";
        return 1;
    }

    sna::Graph<double> graph;
    std::string line;
    int line_number = 0;

    while (std::getline(input, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);

        std::istringstream fields(line);
        int u = 0;
        int v = 0;
        double w = 1.0;
        if (!(fields >> u)) continue;  // blank line
        if (!(fields >> v)) {
            std::cerr << "line " << line_number << ": expected two vertex ids\n";
            return 1;
        }
        if (!(fields >> w)) {
            if (!fields.eof()) {
                std::cerr << "line " << line_number << ": invalid edge weight\n";
                return 1;
            }
            w = 1.0;
        }
        graph.addEdge(u, v, w);
    }

    sna::DetectorParams<double> params;
    if (argc > 2) {
        try {
            params.shift = std::stod(argv[2]);
        } catch (const std::exception&) {
            std::cerr << "invalid shift: " << argv[2] << "\n";
            return 2;
        }
    }

    const std::vector<std::vector<int>> groups = sna::detectSymmetries(graph, params);
    std::cout << "vertices: " << graph.numNodes() << ", edges: " << graph.numEdges()
              << ", shift: " << params.shift << "\n";
    std::cout << "groups: " << groups.size() << "\n";
    for (const std::vector<int>& group : groups) {
        std::cout << "  " << group.size() << " vertices:";
        for (const int vertex : group) std::cout << " " << vertex;
        std::cout << "\n";
    }
    return 0;
}
