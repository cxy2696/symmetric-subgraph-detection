// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Xinyi Cui
//
// Tests for the graph container and the label translation of the detector.

#include <cmath>
#include <tuple>
#include <vector>

#include "sna/sna.hpp"
#include "test_util.hpp"

using namespace sna;

namespace {

void testConstruction() {
    Graph<double> graph;
    CHECK(graph.numNodes() == 0);
    CHECK(graph.labels().empty());

    CHECK(graph.addNode(30));
    CHECK(graph.addNode(10));
    CHECK(!graph.addNode(30));
    CHECK(graph.hasNode(10));
    CHECK(!graph.hasNode(20));
    CHECK(graph.numNodes() == 2);

    const std::vector<int> labels = graph.labels();
    CHECK(labels.size() == 2);
    CHECK(labels[0] == 10);
    CHECK(labels[1] == 30);
}

void testEdges() {
    Graph<double> graph;
    graph.addEdge(5, 7, 2.0);
    graph.addEdge(5, 7, 3.0);  // parallel edges sum
    CHECK(graph.numNodes() == 2);
    CHECK(graph.numEdges() == 1);
    CHECK_NEAR(graph.edgeWeight(5, 7), 5.0, 1e-15);
    CHECK_NEAR(graph.edgeWeight(7, 5), 5.0, 1e-15);
    CHECK(graph.hasEdge(5, 7));
    CHECK(graph.hasEdge(7, 5));
    CHECK(!graph.hasEdge(5, 5));

    graph.addEdge(9, 9);  // self-loop: node created, edge ignored
    CHECK(graph.hasNode(9));
    CHECK(graph.numEdges() == 1);
    CHECK(graph.degree(9) == 0);

    graph.addEdge(5, 9);
    CHECK(graph.degree(5) == 2);
    CHECK(graph.degree(7) == 1);
    CHECK(graph.numEdges() == 2);

    const std::vector<int> neighbors = graph.neighbors(5);
    CHECK(neighbors.size() == 2);
    CHECK(neighbors[0] == 7);
    CHECK(neighbors[1] == 9);

    const auto edges = graph.edges();
    CHECK(edges.size() == 2);
    CHECK(std::get<0>(edges[0]) == 5);
    CHECK(std::get<1>(edges[0]) == 7);
    CHECK_NEAR(std::get<2>(edges[0]), 5.0, 1e-15);
    CHECK(std::get<0>(edges[1]) == 5);
    CHECK(std::get<1>(edges[1]) == 9);
}

void testAdjacencyMatrix() {
    Graph<double> graph;
    graph.addNode(40);
    graph.addEdge(30, 10);
    graph.addEdge(20, 10);

    const std::vector<int> labels = graph.labels();
    CHECK(labels.size() == 4);
    CHECK(labels[0] == 10);
    CHECK(labels[1] == 20);
    CHECK(labels[2] == 30);
    CHECK(labels[3] == 40);

    const Matrix<double> a = graph.adjacencyMatrix();
    CHECK(a.rows() == 4);
    CHECK(a.cols() == 4);
    CHECK(isSymmetric(a));
    CHECK_NEAR(a(0, 1), 1.0, 0.0);  // 10 -- 20
    CHECK_NEAR(a(0, 2), 1.0, 0.0);  // 10 -- 30
    CHECK_NEAR(a(1, 2), 0.0, 0.0);  // 20 not connected to 30
    CHECK_NEAR(a(3, 3), 0.0, 0.0);  // isolated 40
    CHECK(graph.degree(10) == 2);
}

void testDensity() {
    Graph<double> triangle;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) triangle.addEdge(i, j);
    }
    CHECK_NEAR(triangle.density({0, 1, 2}), 1.0, 1e-15);
    CHECK_NEAR(triangle.density({0, 1}), 1.0, 1e-15);
    CHECK_NEAR(triangle.density({0}), 0.0, 0.0);

    Graph<double> path;
    for (int i = 0; i < 3; ++i) path.addEdge(i, i + 1);
    CHECK_NEAR(path.density({0, 1, 2, 3}), 0.5, 1e-15);

    Graph<double> weighted;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) weighted.addEdge(i, j, 2.0);
    }
    CHECK_NEAR(weighted.density({0, 1, 2}), 2.0, 1e-15);

    Graph<double> empty;
    CHECK_NEAR(empty.density({0, 1}), 0.0, 0.0);

    // Duplicate labels are treated as a set.
    CHECK_NEAR(triangle.density({0, 0, 1, 1, 2, 2}), 1.0, 1e-15);
}

void testZeroWeightEdge() {
    Graph<double> graph;
    graph.addEdge(0, 1, 0.0);
    CHECK(graph.hasEdge(0, 1));
    CHECK(graph.numEdges() == 1);
    CHECK_NEAR(graph.edgeWeight(0, 1), 0.0, 0.0);
    CHECK_NEAR(graph.adjacencyMatrix()(0, 1), 0.0, 0.0);
    CHECK_NEAR(graph.density({0, 1}), 0.0, 0.0);
}

void testLabelTranslation() {
    Graph<double> graph;
    graph.addNode(500);
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            graph.addEdge(10 * (i + 1), 10 * (j + 1));
        }
    }
    graph.addEdge(100, 200);

    const std::vector<std::vector<int>> groups = detectSymmetries(graph);
    const std::vector<std::vector<int>> expected = {{10, 20, 30}, {100, 200}};
    CHECK_GROUPS(groups, expected);
    for (const std::vector<int>& group : groups) {
        for (const int label : group) CHECK(graph.hasNode(label));
    }
}

void testConstructorAndEmpty() {
    Graph<double> graph(4);
    CHECK(graph.numNodes() == 4);
    CHECK(graph.numEdges() == 0);
    CHECK(graph.adjacencyMatrix().rows() == 4);
    CHECK(detectSymmetries(graph).empty());

    graph.addEdge(0, 1);
    const std::vector<std::vector<int>> groups = detectSymmetries(graph);
    const std::vector<std::vector<int>> expected = {{0, 1}};
    CHECK_GROUPS(groups, expected);

    Graph<double> empty;
    CHECK(detectSymmetries(empty).empty());
}

}  // namespace

void runGraphTests() {
    std::cout << "test_graph\n";
    testConstruction();
    testEdges();
    testAdjacencyMatrix();
    testDensity();
    testZeroWeightEdge();
    testLabelTranslation();
    testConstructorAndEmpty();
}
