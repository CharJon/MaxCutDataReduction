#include "networkit/graph/Graph.hpp"
#include <gtest/gtest.h>

#include "mcp/solver/enumeration_solver.hpp"

TEST(enumerationSolverNk, simpleNoFixed) {
    NetworKit::Graph g(9, true, false);
    g.addEdge(0, 1, -17);
    g.addEdge(0, 2, -5);
    g.addEdge(1, 2, -3);
    g.addEdge(0, 3, 10);
    g.addEdge(0, 4, 10);
    g.addEdge(1, 5, 10);
    g.addEdge(1, 6, 10);
    g.addEdge(2, 7, 10);
    g.addEdge(2, 8, 10);
    ASSERT_EQ(g.numberOfNodes(), 9);
    ASSERT_EQ(g.numberOfEdges(), 9);

    auto solver = mcp::EnumerationSolverNk(&g);
    solver.solve();
    EXPECT_EQ(solver.bestValue(), 60.0);
}

TEST(enumerationSolverNk, simpleFixed) {
    NetworKit::Graph g(9, true, false);
    g.addEdge(0, 1, -17);
    g.addEdge(0, 2, -5);
    g.addEdge(1, 2, -3);
    g.addEdge(0, 3, 10);
    g.addEdge(0, 4, 10);
    g.addEdge(1, 5, 10);
    g.addEdge(1, 6, 10);
    g.addEdge(2, 7, 10);
    g.addEdge(2, 8, 10);
    ASSERT_EQ(g.numberOfNodes(), 9);
    ASSERT_EQ(g.numberOfEdges(), 9);

    auto solver = mcp::EnumerationSolverNk(&g, {0, 1, 2}, {0, 1, 1});
    solver.solve();
    EXPECT_EQ(solver.bestValue(), 38.0);
    auto bestPart = solver.bestPartitionVec();
    EXPECT_EQ(bestPart.size(), 9);
}
