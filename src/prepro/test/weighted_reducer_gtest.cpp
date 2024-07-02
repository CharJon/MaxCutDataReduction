#include "networkit/graph/Graph.hpp"
#include <gtest/gtest.h>

#include "mcp/prepro/weighted_reducer.hpp"

class WeightedReducerTest : public ::testing::Test {};

TEST_F(WeightedReducerTest, path) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, -1);
    g.addEdge(2, 3, 1);
    ASSERT_EQ(g.numberOfNodes(), 4);

    auto wr = mcp::WeightedReducer(g);
    wr.run();
}

TEST_F(WeightedReducerTest, recursive) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 2);
    g.addEdge(0, 2, 2);
    g.addEdge(1, 2, 2);
    g.addEdge(1, 3, -10);
    g.addEdge(1, 4, 5);
    g.addEdge(2, 3, 9);
    g.addEdge(2, 4, -7);
    g.addEdge(3, 4, 7);
    ASSERT_EQ(g.numberOfNodes(), 5);

    auto wr = mcp::WeightedReducer(g);
    wr.run();
    ASSERT_EQ(g.numberOfNodes(), 0);
}
