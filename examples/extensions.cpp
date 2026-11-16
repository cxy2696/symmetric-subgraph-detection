// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Three uses beyond the default setting:
//   1. weighted cliques, where the repeated eigenvalue is minus the weight;
//   2. several shifts in one call;
//   3. post-processing the detected groups (density filter).

#include <iostream>
#include <vector>

#include "sna/sna.hpp"

namespace {

void print(const char* label, const std::vector<std::vector<int>>& groups) {
    std::cout << label << ":\n";
    for (const std::vector<int>& group : groups) {
        std::cout << "  {";
        for (std::size_t i = 0; i < group.size(); ++i) {
            std::cout << (i == 0 ? "" : ", ") << group[i];
        }
        std::cout << "}\n";
    }
    if (groups.empty()) std::cout << "  (none)\n";
}

}  // namespace

int main() {
    using namespace sna;

    // 1. A K4 with every edge weight 2. Its repeated eigenvalue is -2, so the
    //    default shift -1 is not the right choice.
    Graph<double> weighted;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) weighted.addEdge(i, j, 2.0);
    }
    DetectorParams<double> weighted_params;
    weighted_params.shift = -2.0;
    print("weighted K4, shift -2", detectSymmetries(weighted, weighted_params));

    // 2. A triangle has eigenvalue -1; a 4-cycle has eigenvalue 0. A sweep over
    //    both shifts finds both structures.
    Graph<double> mixed;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) mixed.addEdge(i, j);
    }
    for (int i = 0; i < 4; ++i) mixed.addEdge(10 + i, 10 + (i + 1) % 4);
    const std::vector<std::vector<int>> both =
        detectSymmetries(mixed, std::vector<double>{-1.0, 0.0});
    print("triangle + 4-cycle, shifts -1 and 0", both);

    // 3. The 4-cycle groups are independent sets, with density 0. Keep only
    //    groups that are at least 90% dense.
    print("same groups, density filter 0.9", filterByDensity(mixed, both, 0.9));
    return 0;
}
