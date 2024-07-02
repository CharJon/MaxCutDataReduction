#include "networkit/io/EdgeListReader.hpp"
#include <gtest/gtest.h>

#include "mcp/prepro/spqr.hpp"
#include "mcp/solver/enumeration_solver.hpp"

TEST(SPQR, mo5) {
    NetworKit::EdgeListReader r(' ', 0, "#", true, false);
    auto g = r.read("data/mo.wel");

    auto kernelizer = mcp::SPQR(&g, 5);
    kernelizer.run();
    std::cout << kernelizer.offset() << std::endl;

    auto continuesNodeIds = NetworKit::GraphTools::getContinuousNodeIds(g);
    auto graphCompact = NetworKit::GraphTools::getCompactedGraph(g, continuesNodeIds);
    auto enumSolver = mcp::EnumerationSolverNk(&graphCompact);
    enumSolver.solve();

    ASSERT_EQ(kernelizer.offset() + enumSolver.bestValue(), 63);
}

TEST(SPQR, mo6) {
    NetworKit::EdgeListReader r(' ', 0, "#", true, false);
    auto g = r.read("data/mo.wel");

    auto kernelizer = mcp::SPQR(&g, 6);
    kernelizer.run();

    auto continuesNodeIds = NetworKit::GraphTools::getContinuousNodeIds(g);
    auto graphCompact = NetworKit::GraphTools::getCompactedGraph(g, continuesNodeIds);
    auto enumSolver = mcp::EnumerationSolverNk(&graphCompact);
    enumSolver.solve();

    ASSERT_EQ(kernelizer.offset() + enumSolver.bestValue(), 63);
}

TEST(SPQR, mo7) {
    NetworKit::EdgeListReader r(' ', 0, "#", true, false);
    auto g = r.read("data/mo.wel");

    auto kernelizer = mcp::SPQR(&g, 7);
    kernelizer.run();

    auto continuesNodeIds = NetworKit::GraphTools::getContinuousNodeIds(g);
    auto graphCompact = NetworKit::GraphTools::getCompactedGraph(g, continuesNodeIds);
    auto enumSolver = mcp::EnumerationSolverNk(&graphCompact);
    enumSolver.solve();

    ASSERT_EQ(kernelizer.offset() + enumSolver.bestValue(), 63);
}

TEST(SPQR, mo20) {
    NetworKit::EdgeListReader r(' ', 0, "#", true, false);
    auto g = r.read("data/mo.wel");

    auto kernelizer = mcp::SPQR(&g, 15);
    kernelizer.run();

    auto continuesNodeIds = NetworKit::GraphTools::getContinuousNodeIds(g);
    auto graphCompact = NetworKit::GraphTools::getCompactedGraph(g, continuesNodeIds);
    auto enumSolver = mcp::EnumerationSolverNk(&graphCompact);
    enumSolver.solve();

    ASSERT_EQ(kernelizer.offset() + enumSolver.bestValue(), 63);
}