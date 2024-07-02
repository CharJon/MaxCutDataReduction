#include "networkit/io/EdgeListReader.hpp"
#include <gtest/gtest.h>

#include "mcp/prepro/data_reduction.hpp"

TEST(DataReduction, Triangle1) {
    NetworKit::Graph g(3, true, false);
    g.addEdge(0, 1, -10);
    g.addEdge(0, 2, -10);
    g.addEdge(1, 2, 5);
    ASSERT_EQ(g.numberOfNodes(), 3);
    ASSERT_EQ(g.numberOfEdges(), 3);

    auto kernelizer = mcp::DataReduction(&g, 3);
    kernelizer.run();

    EXPECT_EQ(kernelizer.numNodesResult(), 0);
    EXPECT_EQ(kernelizer.offset(), 0);
}

TEST(DataReduction, Triangle2) {
    NetworKit::Graph g(3, true, false);
    g.addEdge(0, 1, 5);
    g.addEdge(0, 2, 5);
    g.addEdge(1, 2, 20);
    ASSERT_EQ(g.numberOfNodes(), 3);
    ASSERT_EQ(g.numberOfEdges(), 3);

    auto kernelizer = mcp::DataReduction(&g, 3);
    kernelizer.run();

    EXPECT_EQ(kernelizer.numNodesResult(), 0);
    EXPECT_EQ(kernelizer.offset(), 25);
}

TEST(DataReduction, TwoJoinedTrianglesUnitWeights) {
    NetworKit::Graph g(5, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(0, 4, 1);
    g.addEdge(3, 4, 1);
    ASSERT_EQ(g.numberOfNodes(), 5);
    ASSERT_EQ(g.numberOfEdges(), 6);

    auto kernelizer = mcp::DataReduction(&g, 3);
    kernelizer.run();

    EXPECT_EQ(kernelizer.numNodesResult(), 0);
    EXPECT_EQ(kernelizer.offset(), 4);
}
