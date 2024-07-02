#include <random>

#include "networkit/graph/Graph.hpp"
#include <gtest/gtest.h>

#include "mcp/solver/burer_solver.hpp"

TEST(BurerSolver, oneNode) {
    NetworKit::Graph g(2);
    g.addEdge(0, 1);

    auto solver = mcp::BurerMqlibSolver(g, 1);
    auto bestValue = solver.run();

    EXPECT_TRUE(solver.hasRun());
    EXPECT_EQ(bestValue, 1);
}