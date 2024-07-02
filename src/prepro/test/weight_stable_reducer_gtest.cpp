#include "mcp/prepro/weight_stable_reducer.hpp"

#include "networkit/graph/Graph.hpp"
#include <gtest/gtest.h>

TEST(weightStableReducer, triangle) {
    NetworKit::Graph g(3, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(1, 2, 1);
    ASSERT_EQ(g.numberOfNodes(), 3);
    ASSERT_EQ(g.numberOfEdges(), 3);

    auto ef = mcp::WeightStableReducer(g);
    ef.run();
    EXPECT_EQ(g.numberOfNodes(), 0);
    EXPECT_EQ(g.numberOfEdges(), 0);
    EXPECT_EQ(ef.getOffset(), 2);
}

TEST(weightStableReducer, fivePath) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(4, 0, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 5);

    auto ef = mcp::WeightStableReducer(g);
    ef.run();
    EXPECT_EQ(g.numberOfNodes(), 0);
    EXPECT_EQ(g.numberOfEdges(), 0);
    EXPECT_EQ(ef.getOffset(), 4);
}

TEST(weightStableReducer, house) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(4, 0, 1);
    g.addEdge(4, 1, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 6);

    auto ef = mcp::WeightStableReducer(g);
    ef.run();
    EXPECT_EQ(g.numberOfNodes(), 0);
    EXPECT_EQ(g.numberOfEdges(), 0);
    EXPECT_EQ(ef.getOffset(), 5);
}

class WeightStableReducerTest : public mcp::WeightStableReducer {
public:
    WeightStableReducerTest(NetworKit::Graph &g) : mcp::WeightStableReducer(g) {}

    FRIEND_TEST(Kernelizer, cliqueFinderTestClique);
};

TEST(Kernelizer, cliqueFinderTestClique) {
    NetworKit::Graph g(12, true, false);
    for (int i = 0; i < 7; i++) {
        for (int j = i + 1; j < 7; j++) {
            g.addEdge(i, j, 1);
        }
    }
    for (int i = 3; i < 7; i++) {
        for (int j = 7; j < 12; j++) {
            g.addEdge(i, j, 1);
        }
    }
    for (int i = 7; i < 12; i++) {
        for (int j = i + 1; j < 12; j++) {
            g.addEdge(i, j, 1);
        }
    }

    WeightStableReducerTest algo(g);
    auto ret = algo.getOneRemovableClique();
    EXPECT_TRUE(ret.has_value());
}