#include <gtest/gtest.h>
#include "mcp/auxiliary/io.hpp"

#include "networkit/graph/Graph.hpp"

TEST(nkToOgdf, simple) {
    NetworKit::Graph g(4, true, false);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, -1);
    g.addEdge(2, 3, 1);
    ASSERT_EQ(g.numberOfNodes(), 4);

    auto [ogdfGraph, map1, map2] = mcp::nkToOgdf(g);
    ASSERT_EQ(ogdfGraph.numberOfNodes(), 4);
    ASSERT_EQ(ogdfGraph.numberOfEdges(), 3);
}