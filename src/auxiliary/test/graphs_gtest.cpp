#include "networkit/graph/Graph.hpp"
#include <gtest/gtest.h>

#include "mcp/auxiliary/graphs.hpp"

TEST(graphTools, complementGraphLine) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, -1);
    g.addEdge(2, 3, 1);
    ASSERT_EQ(g.numberOfNodes(), 4);

    auto complementGraph = mcp::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 3);
}

TEST(graphTools, complementGraphCycle1) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, -1);
    g.addEdge(2, 4, 1);
    g.addEdge(4, 3, 1);
    g.addEdge(3, 0, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto complementGraph = mcp::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 5);
}

TEST(graphTools, complementGraphCycle2) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 4, 1);
    g.addEdge(4, 1, -1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(3, 0, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto complementGraph = mcp::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 5);
}

TEST(graphTools, complementGraphClique) {
    NetworKit::Graph g(5, true, false);
    for (NetworKit::node u = 0; u < 5; u++) {
        for (NetworKit::node v = u + 1; v < 5; v++) {
            g.addEdge(u, v);
        }
    }
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto complementGraph = mcp::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 0);
}

TEST(graphTools, complementGraphNearClique) {
    NetworKit::Graph g(5, true, false);
    for (NetworKit::node u = 0; u < 5; u++) {
        for (NetworKit::node v = u + 1; v < 5; v++) {
            g.addEdge(u, v);
        }
    }
    g.removeEdge(0, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto complementGraph = mcp::unweightedComplementGraph(g);
    EXPECT_EQ(complementGraph.numberOfNodes(), g.numberOfNodes());
    EXPECT_EQ(complementGraph.numberOfEdges(), 1);
}