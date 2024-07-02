#include <gtest/gtest.h>

#include "mcp/instance/maxcut.hpp"

TEST(MaxCut, createClassEmpty) {
    NetworKit::Graph g(0);
    mcp::MaxCut maxcut(g);
    ASSERT_EQ(maxcut.getNumberOfNodes(), 0);
    ASSERT_EQ(maxcut.getNumberOfEdges(), 0);
    ASSERT_EQ(maxcut.getScalingFactor(), 1);
    ASSERT_EQ(maxcut.getOffset(), 0);
}

TEST(MaxCut, createClassBasic) {
    NetworKit::Graph g(3);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 0);

    mcp::MaxCut maxcut(g);
    ASSERT_EQ(maxcut.getNumberOfNodes(), 3);
    ASSERT_EQ(maxcut.getNumberOfEdges(), 3);
    ASSERT_EQ(maxcut.getScalingFactor(), 1);
    ASSERT_EQ(maxcut.getOffset(), 0);
}

TEST(MaxCut, createClassNonCompact) {
    NetworKit::Graph g(4);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 0);
    g.addEdge(3, 1);
    g.addEdge(3, 2);
    g.removeNode(0);

    mcp::MaxCut maxcut(g);
    ASSERT_EQ(maxcut.getNumberOfNodes(), 3);
    ASSERT_EQ(maxcut.getNumberOfEdges(), 3);
    ASSERT_EQ(maxcut.getScalingFactor(), 1);
    ASSERT_EQ(maxcut.getOffset(), 0);

    ASSERT_EQ(maxcut.getOriginalNode(0), 1);
    ASSERT_EQ(maxcut.getNewNode(1), 0);
}

TEST(MaxCut, solutionValue) {
    NetworKit::Graph g(3, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 2);
    g.addEdge(2, 0, 3);

    mcp::MaxCut maxcut(g, 2, 3);
    ASSERT_EQ(maxcut.getNumberOfNodes(), 3);
    ASSERT_EQ(maxcut.getNumberOfEdges(), 3);
    ASSERT_EQ(maxcut.getScalingFactor(), 2);
    ASSERT_EQ(maxcut.getOffset(), 3);

    EXPECT_EQ(maxcut.getSolutionValue({0, 1, 1}), 4);
    EXPECT_EQ(maxcut.getOriginalSolutionValue({0, 1, 1}), 11);
}
