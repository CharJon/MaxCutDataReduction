#include <random>

#include "networkit/graph/Graph.hpp"
#include <gtest/gtest.h>

#include "mcp/auxiliary/graphs.hpp"
#include "mcp/solver/unit_weight_solver.hpp"

TEST(unit_weight_solver_gtest, tinyTest) {
    NetworKit::Graph g(5);
    g.addEdge(0, 2);
    g.addEdge(0, 3);
    g.addEdge(0, 4);
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(1, 4);
    g.addEdge(0, 1);
    g.addEdge(2, 3);

    NetworKit::Graph gCop = g;

    mcp::UnitWeightSolver solver(g, 1.0, true);
    EXPECT_TRUE(solver.isApplicable());
    solver.solve();

    EXPECT_EQ(solver.bestValue(), 6);

    std::vector<uint8_t> bestSolution = solver.bestPartitionVec();

    int numEdgesInCut = 0;

    for (int i = 0; i < 5; i++) {
        EXPECT_TRUE(bestSolution[i] == 0 || bestSolution[i] == 1);
        for (int j = i + 1; j < 5; j++) {
            if (bestSolution[i] != bestSolution[j]) {
                EXPECT_TRUE(gCop.hasEdge(i, j));
                numEdgesInCut++;
            }
        }
    }

    EXPECT_EQ(numEdgesInCut, 6);
}

TEST(unit_weight_solver_gtest, notApplicable) {
    NetworKit::Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 4);

    mcp::UnitWeightSolver solver(g, 1.0, true);
    EXPECT_FALSE(solver.isApplicable());
    EXPECT_ANY_THROW(solver.solve());
}

TEST(unit_weight_solver_gtest, largeApplicableKhalfhalf) {
    uint n = 100;

    std::mt19937 rng(42);
    std::uniform_int_distribution<uint> diss(0, n);

    NetworKit::Graph g(n);

    int numMissingEdges = 0;
    for (uint i = 0; i < n; i++) {
        for (uint j = i + 1; j < n; j++) {
            if (i < n / 2 && j >= n / 2) {
                g.addEdge(i, j);
            } else {
                if (diss(rng) % 2 == 0)
                    g.addEdge(i, j);
                else
                    numMissingEdges++;
            }
        }
    }

    EXPECT_EQ(numMissingEdges, mcp::unweightedComplementGraph(g).numberOfEdges());

    mcp::UnitWeightSolver solver(g, 1.0, true);
    EXPECT_TRUE(solver.isApplicable());
    solver.solve();

    EXPECT_EQ(solver.bestValue(), (n / 2) * ((n + 1) / 2));

    std::vector<uint8_t> bestSolution = solver.bestPartitionVec();

    int numEdgesInCut = 0;
    for (uint i = 0; i < n; i++) {
        EXPECT_TRUE(bestSolution[i] == 0 || bestSolution[i] == 1);
        for (uint j = i + 1; j < n; j++) {
            if (bestSolution[i] != bestSolution[j]) {
                EXPECT_TRUE(g.hasEdge(i, j));
                numEdgesInCut++;
            }
        }
    }

    EXPECT_EQ(numEdgesInCut, (n / 2) * ((n + 1) / 2));
}

TEST(unit_weight_solver_gtest, largeNonApplicableTest) {
    unsigned int n = 100;

    std::mt19937 rng(42);
    std::uniform_int_distribution<unsigned int> diss(0, n);

    NetworKit::Graph g(n);

    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i + 1; j < n; j++) {
            if (i <= n / 10 || j != i + 1) {
                if (diss(rng) % 4 != 0)
                    g.addEdge(i, j);
            }
        }
    }

    NetworKit::Graph gCop = g;

    mcp::UnitWeightSolver solver(g, 1.0, true);
    EXPECT_FALSE(solver.isApplicable());
    EXPECT_ANY_THROW(solver.solve());

    EXPECT_EQ(g.numberOfNodes(), gCop.numberOfNodes());
    EXPECT_EQ(g.numberOfEdges(), gCop.numberOfEdges());
}