#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>

#include "mcp/auxiliary/block_cut_tree.hpp"


TEST(BlockCutTree, SimpleGraph) {
    NetworKit::Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(0, 3);
    g.addEdge(0, 4);
    g.addEdge(3, 4);
    g.addEdge(1, 2);

    BlockCutTree blockCutTree(g);
    blockCutTree.run();

    ASSERT_EQ(blockCutTree.numberOfComponents(), 2);
    ASSERT_EQ(blockCutTree.getArticulationPoints().size(), 1);
}

TEST(BlockCutTree, SimpleGraphWithBridge) {
    NetworKit::Graph g(6);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);

    g.addEdge(0, 5);

    g.addEdge(5, 4);
    g.addEdge(3, 4);
    g.addEdge(3, 5);

    BlockCutTree blockCutTree(g);
    blockCutTree.run();

    ASSERT_EQ(blockCutTree.numberOfComponents(), 3);
    ASSERT_EQ(blockCutTree.getArticulationPoints().size(), 2);
}

TEST(BlockCutTree, SimpleGraphBigger) {
    NetworKit::Graph g(11);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 4);
    g.addEdge(4, 5);
    g.addEdge(4, 6);
    g.addEdge(4, 7);
    g.addEdge(4, 8);
    g.addEdge(5, 6);
    g.addEdge(7, 9);
    g.addEdge(8, 9);
    g.addEdge(9, 10);

    BlockCutTree blockCutTree(g);
    blockCutTree.run();

    ASSERT_EQ(blockCutTree.numberOfComponents(), 5);
    ASSERT_EQ(blockCutTree.getArticulationPoints().size(), 3);

    ASSERT_FALSE(blockCutTree.isArticulationPoint(0));
    ASSERT_FALSE(blockCutTree.isArticulationPoint(10));
    ASSERT_TRUE(blockCutTree.isArticulationPoint(2));
    ASSERT_TRUE(blockCutTree.isArticulationPoint(4));

    ASSERT_EQ(blockCutTree.getComponentsOfNode(4).size(), 3);
    for (auto &comp: blockCutTree.getComponentsOfNode(4)) {
        EXPECT_THAT(blockCutTree.getArticulationPointsOfComponent(comp), testing::Contains(4));
    }
}
