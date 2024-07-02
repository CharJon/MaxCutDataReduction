#include <chrono>
#include <cstdlib>
#include <iostream>

#include "networkit/auxiliary/Parallelism.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"
// #include "networkit/io/EdgeListWriter.hpp"

#include "mcp/auxiliary/chrono_io.hpp"
#include "mcp/auxiliary/cl_parser.hpp"
#include "mcp/auxiliary/io.hpp"
#include "mcp/auxiliary/statistics.hpp"
#include "mcp/prepro/data_reduction.hpp"
#include "mcp/solver/burer_solver.hpp"
#include "mcp/solver/gurobi_quadratic_solver.hpp"

int main(int argc, char *argv[]) {
    std::cout << "Program started." << std::endl;
    Aux::setNumberOfThreads(1); // Deactivate possible parallelism of NetworKit

    // we are lazy and use the same parser as in the data-reduction
    auto clParser = mcp::CLParser(argc, argv);

    std::cout << "********************************************" << std::endl;
    std::cout << "* BQP solver based on Gurobi               *" << std::endl;
    std::cout << "********************************************" << std::endl << std::endl;

    NetworKit::Graph g = mcp::lazyGraphReader(clParser.getFileName(), clParser.getType());
    if (g.isDirected()) {
        std::cout << "Warning: Graph is directed. Will convert to undirected." << std::endl;
        g = NetworKit::GraphTools::toUndirected(g);
    }

    // self-loops
    std::cout << "Removing " << g.numberOfSelfLoops() << " selfloops." << std::endl;
    g.removeSelfLoops();

    // parallel-loops

    // degree zero nodes
    std::vector<NetworKit::node> degreeZeroNodes{};
    for (auto v : g.nodeRange()) {
        if (g.degree(v) == 0)
            degreeZeroNodes.push_back(v);
    }
    std::cout << "Removing " << degreeZeroNodes.size() << " degree zero nodes." << std::endl;
    for (auto v : degreeZeroNodes) {
        g.removeNode(v);
    }

    if (clParser.getSeed() == 0) {
        auto a = NetworKit::GraphTools::getContinuousNodeIds(g);
        g = NetworKit::GraphTools::getCompactedGraph(g, a);
    } else if (clParser.getSeed() >= 1) {
        Aux::Random::setSeed(clParser.getSeed(), false);
        auto a = NetworKit::GraphTools::getRandomContinuousNodeIds(g);
        g = NetworKit::GraphTools::getCompactedGraph(g, a);
    }

    if (!g.isWeighted()) {
        std::cout << "Warning: Graph is unweighted. Will convert to weighted." << std::endl;
        g = NetworKit::GraphTools::toWeighted(g);
    }

    if (clParser.getScalingFactor().has_value()) {
        mcp::scaleEdgeWeights(g, clParser.getScalingFactor().value());
    }

    // print stats
    std::cout << "---------- Graph statistics -----------------" << std::endl;
    auto graphStats = mcp::graphMetricsJson(g, true);
    std::cout << graphStats.dump(4) << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    mcp::simpleMcIntWriter(g, clParser.getFileName().append("_clean.mc"));

    auto outputJson = nlohmann::json();
    outputJson["original_graph_stats"] = graphStats;

    if (clParser.getSolve()) {
        if (clParser.getConnectivity() >= 0) {

            std::cout << "Apply reduction with connectivity " << clParser.getConnectivity() << std::endl;
            auto numNodesInput = g.numberOfNodes();
            auto numEdgesInput = g.numberOfEdges();

            auto reducer = mcp::DataReduction(&g, clParser.getConnectivity(), clParser.getTriangleFlag(),
                                              clParser.getCliqueLevel(), clParser.getStableReduceFlag());
            reducer.run();
            double runtimeReduction = reducer.getRuntimeInSeconds();
            double solutionValue = reducer.offset();

            outputJson["reduction_stats"] = reducer.getStatsJson();
            double totalRuntimeGurobi = 0.;
            bool isOptimal = false;
            int i = 0;
            auto solvingStatsArray = nlohmann::json::array();
            std::cout << "Solving " << reducer.finishedGraphs().size() << " independent graphs." << std::endl;
            for (auto &currentGraph : reducer.finishedGraphs()) {
                nlohmann::json currentGraphJson = {{"graph_id", i}};

                auto a = NetworKit::GraphTools::getContinuousNodeIds(currentGraph);
                auto compactGraph = NetworKit::GraphTools::getCompactedGraph(currentGraph, a);
                currentGraphJson["graph_stats"] = mcp::graphMetricsJson(compactGraph, true);
                // write graph to file
                // NetworKit::EdgeListWriter writer(' ', 0, false);
                // writer.write(compactGraph, "graph_" + std::to_string(i) + ".wel");

                std::cout << "Solving number " << i << " with " << compactGraph.numberOfNodes() << " nodes."
                          << std::endl;
                GurobiQuadraticSolver solver(compactGraph, clParser.getTimeLimit(), clParser.getSeed());
                solutionValue += solver.run();
                currentGraphJson["solver_stats"] = solver.getStats();
                totalRuntimeGurobi += solver.getRuntimeSeconds();
                isOptimal = i == 0 ? solver.isOptimal() : isOptimal && solver.isOptimal();

                solvingStatsArray.emplace_back(currentGraphJson);
                i += 1;
            }
            outputJson["solving_stats"] = solvingStatsArray;

            outputJson["total_runtime_seconds"] = totalRuntimeGurobi + runtimeReduction;
            outputJson["reduction_runtime_seconds"] = runtimeReduction;
            outputJson["gurobi_runtime_seconds"] = totalRuntimeGurobi;
            outputJson["best_solution_value"] = solutionValue;
            outputJson["reduction_offset"] = reducer.offset();
            outputJson["optimal_solution_found"] = isOptimal;
            outputJson["num_independent_graphs"] = reducer.finishedGraphs().size();
            outputJson["seed"] = clParser.getSeed();
            outputJson["node_efficiency"] =
                1 - (static_cast<double>(reducer.numNodesResult()) / static_cast<double>(numNodesInput));
            outputJson["edge_efficiency"] =
                1 - (static_cast<double>(reducer.numEdgesResult()) / static_cast<double>(numEdgesInput));
        } else {
            std::cout << "Connectivity is " << clParser.getConnectivity() << ". Will solve original graph via gurobi."
                      << std::endl;
            GurobiQuadraticSolver solver(g, clParser.getTimeLimit(), clParser.getSeed());
            auto solutionValue = solver.run();
            outputJson["best_solution_value"] = solutionValue;
            outputJson["optimal_solution_found"] = solver.isOptimal();
            outputJson["total_runtime_seconds"] = solver.getRuntimeSeconds();
            outputJson["seed"] = clParser.getSeed();
            outputJson["gurobi_runtime_seconds"] = solver.getRuntimeSeconds();
        }
    } else if (clParser.getHeuristic()) {
        if (clParser.getConnectivity() >= 0) {
            std::cout << "Apply reduction with connectivity " << clParser.getConnectivity() << std::endl;
            auto reducer = mcp::DataReduction(&g, clParser.getConnectivity(), clParser.getTriangleFlag(),
                                              clParser.getCliqueLevel(), clParser.getStableReduceFlag());
            reducer.run();
            auto reducerRuntimeSeconds = reducer.getRuntimeInSeconds();
            double solverRuntimeSeconds = 0.;
            double bestValue = reducer.offset();
            double timeToBestSeconds = 0;
            outputJson["reduction_stats"] = reducer.getStatsJson();

            auto solvingStatsArray = nlohmann::json::array();
            std::cout << "Solving " << reducer.finishedGraphs().size() << " independent graphs." << std::endl;

            if (reducer.finishedGraphs().size() > 0) {
                auto solver = mcp::BurerMqlibSolver(g, clParser.getTimeLimit() - reducerRuntimeSeconds);
                // std::cout << "Best value on orig:   " << solver.run() << ". Took " <<
                // std::to_string(solver.elapsed().count())
                //           << "ms." << std::endl;

                solver = mcp::BurerMqlibSolver(reducer.mergeReducedGraphs(), clParser.getTimeLimit());
                bestValue += solver.run();
                std::cout << "Best value on merged: " << bestValue << ". Took " << solver.elapsed() << std::endl;
                solverRuntimeSeconds =
                    std::chrono::duration_cast<std::chrono::duration<double>>(solver.elapsed()).count();
                timeToBestSeconds = solver.timeToBestSeconds();
                nlohmann::json solvingStats = {};
                solvingStats["best_solution_values"] = solver.pastSolutionValues();
                solvingStats["best_solution_times"] = solver.pastSolutionTimes();
                outputJson["solving_stats"] = solvingStats;
            }

            outputJson["total_runtime_seconds"] = reducerRuntimeSeconds + solverRuntimeSeconds;
            outputJson["reduction_runtime_seconds"] = reducerRuntimeSeconds;
            outputJson["heuristic_runtime_seconds"] = solverRuntimeSeconds;
            outputJson["best_solution_value"] = bestValue;
            outputJson["time_to_best_seconds"] = timeToBestSeconds;
            outputJson["reduction_offset"] = reducer.offset();
            outputJson["optimal_solution_found"] = false;
            outputJson["num_independent_graphs"] = reducer.finishedGraphs().size();
            outputJson["seed"] = clParser.getSeed();
            outputJson["node_efficiency"] =
                1 - (static_cast<double>(reducer.numNodesResult()) / static_cast<double>(g.numberOfNodes()));
            outputJson["edge_efficiency"] =
                1 - (static_cast<double>(reducer.numEdgesResult()) / static_cast<double>(g.numberOfEdges()));
        } else {
            auto solver = mcp::BurerMqlibSolver(g, clParser.getTimeLimit());
            auto solutionValue = solver.run();
            outputJson["best_solution_value"] = solutionValue;
            outputJson["optimal_solution_found"] = false;
            outputJson["seed"] = clParser.getSeed();
            outputJson["heuristic_runtime_seconds"] =
                std::chrono::duration_cast<std::chrono::duration<double>>(solver.elapsed()).count();
            outputJson["total_runtime_seconds"] = outputJson["heuristic_runtime_seconds"];
            nlohmann::json solvingStats = {};
            solvingStats["best_solution_values"] = solver.pastSolutionValues();
            solvingStats["best_solution_times"] = solver.pastSolutionTimes();
            outputJson["solving_stats"] = solvingStats;
        }
    } else {
        std::cout << "Apply reduction with connectivity " << clParser.getConnectivity() << std::endl;
        auto reducer = mcp::DataReduction(&g, clParser.getConnectivity(), clParser.getTriangleFlag(),
                                          clParser.getCliqueLevel(), clParser.getStableReduceFlag());
        reducer.run();
        outputJson["reduction_stats"] = reducer.getStatsJson();

        outputJson["reduction_runtime_seconds"] = reducer.getRuntimeInSeconds();
        outputJson["reduction_offset"] = reducer.offset();
        outputJson["optimal_solution_found"] = false;
        outputJson["num_independent_graphs"] = reducer.finishedGraphs().size();
        outputJson["seed"] = clParser.getSeed();
        outputJson["node_efficiency"] =
            1 - (static_cast<double>(reducer.numNodesResult()) / static_cast<double>(g.numberOfNodes()));
        outputJson["edge_efficiency"] =
            1 - (static_cast<double>(reducer.numEdgesResult()) / static_cast<double>(g.numberOfEdges()));
    }

    auto combinedStats = nlohmann::json();
    combinedStats["graph_statistics"] = graphStats;
    combinedStats["solver_statistics"] = outputJson;
    std::cout << std::setw(4) << outputJson << std::endl;
    if (clParser.getStatsFile().has_value()) {
        std::ofstream o(clParser.getStatsFile().value());
        o << std::setw(4) << combinedStats << std::endl;
    }

    return EXIT_SUCCESS;
}
