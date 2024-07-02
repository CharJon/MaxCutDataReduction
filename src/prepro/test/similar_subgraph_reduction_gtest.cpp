#include "networkit/graph/Graph.hpp"
#include <gtest/gtest.h>

#include "mcp/prepro/similar_subgraph_reduction.hpp"

TEST(similarSubgraphReduction, singleTriangle1) {
    NetworKit::Graph g(3, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 3);
    g.addEdge(1, 2, 1);
    ASSERT_EQ(g.numberOfNodes(), 3);
    ASSERT_EQ(g.numberOfEdges(), 3);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    EXPECT_EQ(ssr.haveSameWeightedNeighborhood(0, 2, true), 1);
}

TEST(similarSubgraphReduction, singleTriangle2) {
    NetworKit::Graph g(3, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 3);
    g.addEdge(1, 2, 3);
    ASSERT_EQ(g.numberOfNodes(), 3);
    ASSERT_EQ(g.numberOfEdges(), 3);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    EXPECT_EQ(ssr.haveSameWeightedNeighborhood(0, 1, true), 1);
}

TEST(similarSubgraphReduction, singleTriangleMergeNeg) {
    NetworKit::Graph g(3, true, false);
    g.addEdge(0, 1, -1);
    g.addEdge(0, 2, 3);
    g.addEdge(1, 2, 3);
    ASSERT_EQ(g.numberOfNodes(), 3);
    ASSERT_EQ(g.numberOfEdges(), 3);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    EXPECT_EQ(ssr.haveSameWeightedNeighborhood(0, 1, true), 1);
    auto edge = ssr.connectedSameNeighborhood();
    ASSERT_TRUE(edge.has_value());
    EXPECT_TRUE(edge.has_value() && (edge.value()[0] == 0 || edge.value()[0] == 1));
    EXPECT_TRUE(edge.has_value() && (edge.value()[1] == 0 || edge.value()[1] == 1));
}

TEST(similarSubgraphReduction, singleTriangleMergePos) {
    NetworKit::Graph g(3, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 3);
    g.addEdge(1, 2, -3);
    ASSERT_EQ(g.numberOfNodes(), 3);
    ASSERT_EQ(g.numberOfEdges(), 3);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    EXPECT_EQ(ssr.haveSameWeightedNeighborhood(0, 1, true), -1);
    auto edge = ssr.connectedSameNeighborhood();
    ASSERT_TRUE(edge.has_value());
    EXPECT_TRUE(edge.value()[0] == 0 || edge.value()[0] == 1);
    EXPECT_TRUE(edge.value()[1] == 0 || edge.value()[1] == 1);
}

TEST(similarSubgraphReduction, doubleTriangle) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(1, 3, 1);
    ASSERT_EQ(g.numberOfNodes(), 4);
    ASSERT_EQ(g.numberOfEdges(), 5);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    EXPECT_TRUE(ssr.haveSameNeighborhood(2, 3));
    EXPECT_TRUE(ssr.haveSameWeightedNeighborhood(2, 3, false));
    EXPECT_FALSE(ssr.haveSameWeightedNeighborhood(0, 1, false));
    EXPECT_TRUE(ssr.haveSameWeightedNeighborhood(0, 1, true));
}

TEST(similarSubgraphReduction, tripleTriangle) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(0, 4, 1);
    g.addEdge(1, 4, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 7);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    // auto nodesWithSameNeighborhood = ssr.findNodesWithSameNeighborhood();
    // ASSERT_EQ(nodesWithSameNeighborhood.size(), 1);
    // EXPECT_EQ(nodesWithSameNeighborhood[0].size(), 3);
}

TEST(similarSubgraphReduction, fiveClique) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(0, 4, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(1, 4, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(2, 4, 1);
    g.addEdge(3, 4, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 10);
    ASSERT_EQ(g.totalEdgeWeight(), 10);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    for (NetworKit::node i = 0; i < g.numberOfNodes(); i++) {
        for (NetworKit::node j = 0; j < g.numberOfNodes(); j++) {
            if (i != j) {
                EXPECT_EQ(ssr.haveSameNeighborhood(i, j), i == j);
                EXPECT_TRUE(ssr.haveSameWeightedNeighborhood(i, j, true));
                EXPECT_EQ(ssr.areTwins(i, j), 1);
            }
        }
    }

    ssr.mergeTwins(0, 1);
    EXPECT_EQ(g.numberOfNodes(), 4);
    EXPECT_EQ(g.numberOfEdges(), 6);
    EXPECT_EQ(g.totalEdgeWeight(), 9);
}

TEST(similarSubgraphReduction, fiveCliqueRun) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(0, 4, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(1, 4, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(2, 4, 1);
    g.addEdge(3, 4, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 10);
    ASSERT_EQ(g.totalEdgeWeight(), 10);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    ssr.run();
    // EXPECT_EQ(g.numberOfNodes(), 4);
    // EXPECT_EQ(g.numberOfEdges(), 6);
    // EXPECT_EQ(g.totalEdgeWeight(), 9);
}

TEST(similarSubgraphReduction, negativeAlphaNotConnected) {
    NetworKit::Graph g(6, true, false);
    g.addEdge(0, 2, 2);
    g.addEdge(0, 3, 3);
    g.addEdge(1, 2, -2);
    g.addEdge(1, 3, -3);
    g.addEdge(2, 4, 1);
    g.addEdge(3, 5, 1);

    ASSERT_EQ(g.numberOfNodes(), 6);
    ASSERT_EQ(g.numberOfEdges(), 6);
    ASSERT_EQ(g.totalEdgeWeight(), 2);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    ssr.run();
    EXPECT_EQ(ssr.getOffset(), -5);
    EXPECT_EQ(g.numberOfNodes(), 5);
    EXPECT_EQ(g.numberOfEdges(), 4);
    EXPECT_EQ(g.totalEdgeWeight(), 12);
}

TEST(similarSubgraphReduction, positiveAlphaConntected) {
    NetworKit::Graph g(6, true, false);
    g.addEdge(0, 2, 2);
    g.addEdge(0, 3, 3);
    g.addEdge(1, 2, 2);
    g.addEdge(1, 3, 3);
    g.addEdge(0, 1, -1);
    g.addEdge(2, 4, 1);
    g.addEdge(3, 5, 1);

    ASSERT_EQ(g.numberOfNodes(), 6);
    ASSERT_EQ(g.numberOfEdges(), 7);
    ASSERT_EQ(g.totalEdgeWeight(), 11);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    ssr.run();
    EXPECT_EQ(ssr.getOffset(), 0);
    EXPECT_EQ(g.numberOfNodes(), 5);
    EXPECT_EQ(g.numberOfEdges(), 4);
    EXPECT_EQ(g.totalEdgeWeight(), 12);
}

TEST(similarSubgraphReduction, negativeAlphaConntected) {
    NetworKit::Graph g(6, true, false);
    g.addEdge(0, 2, 2);
    g.addEdge(0, 3, 3);
    g.addEdge(1, 2, -2);
    g.addEdge(1, 3, -3);
    g.addEdge(0, 1, 4);
    g.addEdge(2, 4, 1);
    g.addEdge(3, 5, 1);

    ASSERT_EQ(g.numberOfNodes(), 6);
    ASSERT_EQ(g.numberOfEdges(), 7);
    ASSERT_EQ(g.totalEdgeWeight(), 6);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    ssr.run();
    EXPECT_EQ(ssr.getOffset(), -1);
    EXPECT_EQ(g.numberOfNodes(), 5);
    EXPECT_EQ(g.numberOfEdges(), 4);
    EXPECT_EQ(g.totalEdgeWeight(), 12);
}

TEST(similarSubgraphReduction, negativeAlpha2Conntected) {
    NetworKit::Graph g(6, true, false);
    g.addEdge(0, 2, 2);
    g.addEdge(0, 3, 3);
    g.addEdge(1, 2, -4);
    g.addEdge(1, 3, -6);
    g.addEdge(0, 1, 4);
    g.addEdge(2, 4, 1);
    g.addEdge(3, 5, 1);

    ASSERT_EQ(g.numberOfNodes(), 6);
    ASSERT_EQ(g.numberOfEdges(), 7);
    ASSERT_EQ(g.totalEdgeWeight(), 1);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    ssr.run();
    EXPECT_EQ(ssr.getOffset(), -6);
    EXPECT_EQ(g.numberOfNodes(), 5);
    EXPECT_EQ(g.numberOfEdges(), 4);
    EXPECT_EQ(g.totalEdgeWeight(), 17);
}

TEST(similarSubgraphReduction, fiveCliqueDiffWeights) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 0.5);
    g.addEdge(0, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(0, 4, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(1, 4, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(2, 4, 1);
    g.addEdge(3, 4, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 10);
    ASSERT_EQ(g.totalEdgeWeight(), 9.5);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    EXPECT_FALSE(ssr.haveSameNeighborhood(0, 1));
    EXPECT_TRUE(ssr.haveSameWeightedNeighborhood(0, 1, true));
    EXPECT_EQ(ssr.areTwins(0, 1), 1);

    ssr.mergeTwins(0, 1);
    EXPECT_EQ(g.numberOfNodes(), 4);
    EXPECT_EQ(g.numberOfEdges(), 6);
    EXPECT_EQ(g.totalEdgeWeight(), 9);
}

TEST(similarSubgraphReduction, fiveCliqueDiffWeightsRun) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 0.5);
    g.addEdge(0, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(0, 4, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(1, 4, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(2, 4, 1);
    g.addEdge(3, 4, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 10);
    ASSERT_EQ(g.totalEdgeWeight(), 9.5);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    ssr.run();
    // EXPECT_EQ(g.numberOfNodes(), 4);
    // EXPECT_EQ(g.numberOfEdges(), 6);
    // EXPECT_EQ(g.totalEdgeWeight(), 9);
}

TEST(similarSubgraphReduction, big) {
    NetworKit::Graph g(13, true, false);
    // 0 - 4 are random nodes, to make 5 - 9 have different neighborhoods
    g.addEdge(0, 5, 11);
    g.addEdge(1, 6, 12);
    g.addEdge(2, 7, 13);
    g.addEdge(3, 8, 14);
    g.addEdge(4, 9, 15);

    g.addEdge(5, 10, 1);
    g.addEdge(5, 11, 2);
    g.addEdge(5, 12, 2);

    g.addEdge(6, 10, 1);
    g.addEdge(6, 11, 2);
    g.addEdge(6, 12, 2);

    g.addEdge(7, 10, 1);
    g.addEdge(7, 11, 2);
    g.addEdge(7, 12, 2);

    g.addEdge(8, 10, 1);
    g.addEdge(8, 11, 2);
    g.addEdge(8, 12, 2);

    g.addEdge(9, 10, 1);
    g.addEdge(9, 11, 2);
    g.addEdge(9, 12, 2);

    ASSERT_EQ(g.numberOfNodes(), 13);
    ASSERT_EQ(g.numberOfEdges(), 20);
    ASSERT_EQ(g.totalEdgeWeight(), 90);

    auto ssr = mcp::SimilarSubgraphReducer(g);
    ssr.run();
    EXPECT_EQ(g.numberOfNodes(), 11);
    EXPECT_EQ(g.numberOfEdges(), 10);
    EXPECT_EQ(g.totalEdgeWeight(), 90);
}
