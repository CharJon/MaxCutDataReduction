#ifndef MCP_DATA_REDUCTION_HPP
#define MCP_DATA_REDUCTION_HPP

#include "similar_subgraph_reduction.hpp"

#include "networkit/auxiliary/Timer.hpp"
#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "nlohmann/json.hpp"

#include "mcp/auxiliary/graphs.hpp"
#include "mcp/prepro/base.hpp"
#include "mcp/prepro/spqr.hpp"
#include "mcp/prepro/weight_stable_reducer.hpp"
#include "mcp/prepro/weighted_reducer.hpp"
#include "mcp/solver/enumeration_solver.hpp"
#include "mcp/solver/gurobi_quadratic_solver.hpp"
#include "mcp/solver/unit_weight_solver.hpp"

namespace mcp {

class DataReduction {
    struct QueueElement {
        NetworKit::Graph graph;
        int guaranteedConnectivity;
    };

public:
    explicit DataReduction(NetworKit::Graph const *graph, int maxConnectivityToUse = 0, int numTriangles = 7,
                           int cliqueLevel = 2, bool useStableReduce = true, unsigned int enumThreshold = 20,
                           bool requireFullSolve = false, int guaranteedConnectivity = 0)
        : graph_(graph),
          maxConnectivityToUse_(maxConnectivityToUse),
          numTrianglesToUse_(numTriangles),
          cliqueLevel_(cliqueLevel),
          useStableReduce_(useStableReduce),
          enumThreshold_(enumThreshold),
          requireFullSolve_(requireFullSolve) {
        auto copyGraph = NetworKit::Graph(*graph);
        removeZeroWeightAndZeroDegree(copyGraph);
        openGraphs_.emplace_back(std::move(copyGraph), guaranteedConnectivity);
        initStats();
    }

    void initStats() {
        stats_["numDegreeZeroContractions"] = 0;
        stats_["numDegreeOneContractions"] = 0;
        stats_["numDegreeTwoContractions"] = 0;
        stats_["numDegreeThreeContractions"] = 0;
        stats_["numHeavyEdgeContractions"] = 0;
        stats_["numTriangleOneContractions"] = 0;
        stats_["numTriangleTwoContractions"] = 0;
        stats_["numTriangleThreeContractions"] = 0;
        stats_["numSpqrLeaveReductions"] = 0;
        stats_["numSpqrNodesRemoved"] = 0;
        stats_["numTriangleThreeContractions"] = 0;
        stats_["numThreePathsContractions"] = 0;
        stats_["numCliqueRemoval"] = 0;
        stats_["numCliqueNodesRemoved"] = 0;
        stats_["numNearCliqueRemoval"] = 0;
        stats_["numNearCliqueNodesRemoved"] = 0;
        stats_["totalRuntimeSeconds"] = 0.0;
    }

    void run() {
        const auto start = std::chrono::high_resolution_clock::now();
        hasRun_ = false;

        while (!openGraphs_.empty()) {
            auto next = std::move(openGraphs_.back());
            openGraphs_.pop_back();
            queuePopReport(next);
            assert(next.graph.checkConsistency());

            if ((numWeightEdges(next.graph, 1) == next.graph.numberOfEdges())
                && tryUnitWeightSolver(next.graph, false)) {
                continue;
            } else

                // graph is small and therefore easy
                if ((next.graph.numberOfNodes() <= enumThreshold_)) {
                    auto res = tryEnumerationSolve(next.graph);
                    assert(res);
                    // offset_ += gurobiSolve(next.graph);
                }
                // graph might be disconnected
                else if ((next.guaranteedConnectivity == 0) && (maxConnectivityToUse_ >= 0)) {
                    if (!tryConnectedSplit(next)) {
                        openGraphs_.emplace_back(std::move(next.graph), 1);
                    }
                }
                // graph might not be biconnected
                else if (next.guaranteedConnectivity == 1 && (maxConnectivityToUse_ >= 1)) {
                    if (!tryBiconSplit(next)) {
                        openGraphs_.emplace_back(std::move(next.graph), 2);
                    }
                }
                // graph is at least biconnected or maxConnectivityToUse_ is zero
                else {
                    if (!solveIndependentComp(next)) {
                        // if we reach this point, and we could not solve the graph, we have to give up
                        std::cout << "Remaining graph is too big to solve with " << next.graph.numberOfNodes() << "."
                                  << std::endl;
                        finishedGraphs_.push_back(std::move(next.graph));
                        if (requireFullSolve_) {
                            // not being able to solve the graph, means we cannot fully solve the instance
                            return;
                        }
                    }
                }
        }
        solvingTime_ =
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
        hasRun_ = true;
    }

    void queuePopReport(const QueueElement &next) const {
        const unsigned int n = next.graph.numberOfNodes();
        const unsigned int m = next.graph.numberOfEdges();
        std::cout << "Popped graph with " << n << " nodes and " << m << " edges"
                  << " (" << (n * n - n) / 2 - m << " missing edges) "
                  << "and guaranteed connectivity of " << next.guaranteedConnectivity << "." << '\n'
                  << "(#open, #unsolved): (" << openGraphs_.size() << ", " << finishedGraphs_.size() << ")"
                  << std::endl;
    }

    /**
     *
     * @param next
     * @return true if the graph was solved or split, false else
     */
    bool splitOrEnumerate(QueueElement &next) {
        if (next.graph.numberOfNodes() == 0) {
            return true;
        }

        if (next.graph.numberOfNodes() <= enumThreshold_) {
            if (numWeightEdges(next.graph, 1.) == next.graph.numberOfEdges())
                if (tryUnitWeightSolver(next.graph, true))
                    return true;
            tryEnumerationSolve(next.graph);
            return true;
        }

        if ((maxConnectivityToUse_ >= 0) && tryConnectedSplit(next)) {
            std::cout << "Connected split successful" << std::endl;
            return true;
        }

        if ((maxConnectivityToUse_ >= 1) && tryBiconSplit(next)) {
            std::cout << "Bicon split successful" << std::endl;
            return true;
        }

        return false;
    }

    bool solveIndependentComp(QueueElement &next) {
        std::vector<NetworKit::node> activeNodes;

        bool graphDidShrink = true;
        unsigned int maxNumIterations = 10;
        for (unsigned int i = 0; (i < maxNumIterations) && graphDidShrink; i++) {
            std::cout << "Data-reducer starting round " << i << std::endl;

            // if graph has unit weights
            if (useStableReduce_ && (numWeightEdges(next.graph, 1.) == next.graph.numberOfEdges())) {
                auto wsr = WeightStableReducer(next.graph, cliqueLevel_ > 1, cliqueLevel_ > 0);
                wsr.run();
                offset_ += wsr.getOffset();
                wsr.sumStats(stats_);
            }

            auto wr = i > 0 ? WeightedReducer(next.graph, activeNodes, (maxConnectivityToUse_ >= 3), numTrianglesToUse_)
                            : WeightedReducer(next.graph, (maxConnectivityToUse_ >= 3), numTrianglesToUse_);
            wr.run();
            offset_ += wr.getOffset();
            wr.sumStats(stats_);

            // ogdf crashes if the graph is not biconnected
            if (tryBiconSplit(next)) {
                return true;
            }

            if ((i > 0) && (next.graph.numberOfNodes() > 1000000)) {
                graphDidShrink = false;
            } else {
                const auto numNodesBefore = next.graph.numberOfNodes();
                if (next.graph.numberOfNodes() > 0 && (maxConnectivityToUse_ >= 2)) {
                    auto spqr = SPQR(&next.graph, enumThreshold_);
                    spqr.run();
                    offset_ += spqr.offset();
                    stats_["numSpqrLeaveReductions"] =
                        static_cast<unsigned int>(stats_["numSpqrLeaveReductions"]) + spqr.numLeafContractions();
                    stats_["numSpqrNodesRemoved"] =
                        static_cast<unsigned int>(stats_["numSpqrNodesRemoved"]) + spqr.numRemovedNodes();
                    activeNodes = spqr.activatedNodes();
                }

                if (i == 0) {
                    auto sr = SimilarSubgraphReducer(next.graph);
                    sr.run();
                    offset_ += sr.getOffset();
                }

                graphDidShrink = (numNodesBefore != next.graph.numberOfNodes());
            }
        }

        return splitOrEnumerate(next);
    }

    bool tryConnectedSplit(const QueueElement &next) {
        NetworKit::ConnectedComponents connectedComponents(next.graph);
        connectedComponents.run();
        if (connectedComponents.numberOfComponents() == 1) {
            return false;
        }

        int numSmallComponents = 0;
        for (const auto &c : connectedComponents.getComponents()) {
            if (c.size() <= 1) {
                numSmallComponents++;
            } else if (c.size() == 2) {
                // handle bridges right away
                auto u = c[0];
                auto v = c[1];
                assert(next.graph.hasEdge(u, v));
                auto weight = next.graph.weight(u, v);
                offset_ += std::max(0., weight);
                numSmallComponents++;
            } else {
                auto inducedSubgraph = NetworKit::GraphTools::subgraphFromNodes(next.graph, c.begin(), c.end(), true);
                openGraphs_.emplace_back(std::move(inducedSubgraph), 1);
            }
        }
        std::cout << "Connected split successful with " << numSmallComponents << " comps <= 2" << std::endl;
        return true;
    }

    bool tryBiconSplit(const QueueElement &next) {
        NetworKit::BiconnectedComponents biconnectedComponents(next.graph);
        biconnectedComponents.run();
        if (biconnectedComponents.numberOfComponents() == 1) {
            return false;
        }
        int numBridges = 0;
        size_t maxCompSize = 0;
        size_t totalSize = 0;
        for (const auto &c : biconnectedComponents.getComponents()) {
            maxCompSize = std::max(maxCompSize, c.size());
            totalSize += c.size() - 1;
            if (c.size() == 2) {
                // handle bridges right away
                auto u = c[0];
                auto v = c[1];
                assert(next.graph.hasEdge(u, v));
                auto weight = next.graph.weight(u, v);
                offset_ += std::max(0., weight);
                numBridges++;
            } else {
                auto inducedSubgraph = NetworKit::GraphTools::subgraphFromNodes(next.graph, c.begin(), c.end(), true);
                openGraphs_.emplace_back(std::move(inducedSubgraph), 2);
            }
        }
        std::cout << "Bicon split successful with " << numBridges << " bridges and "
                  << biconnectedComponents.getComponents().size() << " components. Max size = " << maxCompSize
                  << " total size = " << totalSize << std::endl;
        return true;
    }

    int numNodesResult() const {
        assert(hasRun_);
        if (finishedGraphs_.empty())
            return 0;
        int sumOfNodeNumbers =
            std::accumulate(finishedGraphs_.begin(), finishedGraphs_.end(), 0,
                            [](int cur, NetworKit::Graph const &g) { return cur + g.numberOfNodes(); });
        return sumOfNodeNumbers + 1 - static_cast<int>(finishedGraphs_.size());
    }

    int numEdgesResult() const {
        assert(hasRun_);
        int sumOfEdgeNumbers =
            std::accumulate(finishedGraphs_.begin(), finishedGraphs_.end(), 0,
                            [](int cur, NetworKit::Graph const &g) { return cur + g.numberOfEdges(); });
        return sumOfEdgeNumbers;
    }

    NetworKit::Graph mergeReducedGraphs() const {
        assert(hasRun_);
        NetworKit::Graph mergedGraph(0, true, false);
        for (auto &g : finishedGraphs_) {
            NetworKit::GraphTools::append(mergedGraph, g);
        }
        return mergedGraph;
    }

    unsigned int numFinishedGraphs() const {
        assert(hasRun_);
        return finishedGraphs_.size();
    }

    const auto &finishedGraphs() const {
        assert(hasRun_);
        return finishedGraphs_;
    }

    double offset() const {
        assert(hasRun_);
        return offset_;
    }

    double getRuntimeInSeconds() const {
        assert(hasRun_);
        return std::chrono::duration_cast<std::chrono::duration<double>>(solvingTime_).count();
    }

    nlohmann::json getStatsJson() const {
        assert(hasRun_);
        nlohmann::ordered_json stats;
        stats["numNodesOriginal"] = graph_->numberOfNodes();
        stats["numEdgesOriginal"] = graph_->numberOfEdges();
        stats["numNodesResult"] = numNodesResult();
        stats["numEdgesResult"] = numEdgesResult();
        stats["numFinishedGraphs"] = numFinishedGraphs();
        stats["offset"] = offset();
        stats["runtimeInSeconds"] = getRuntimeInSeconds();
        stats["reducerStats"] = stats_;
        return stats;
    }

private:
    bool tryUnitWeightSolver(NetworKit::Graph const &graph, bool fullTest = false) {
        if (graph.numberOfNodes() != graph.upperNodeIdBound()) {
            auto nodeIdMap = NetworKit::GraphTools::getContinuousNodeIds(graph);
            auto compactG = NetworKit::GraphTools::getCompactedGraph(graph, nodeIdMap);

            auto uws = UnitWeightSolver(compactG, 1., fullTest);
            if (uws.isApplicable()) {
                uws.solve();
                offset_ += uws.bestValue();
                return true;
            }
        } else {
            auto uws = UnitWeightSolver(graph, 1., fullTest);
            if (uws.isApplicable()) {
                uws.solve();
                offset_ += uws.bestValue();
                return true;
            }
        }

        return false;
    }

    bool tryEnumerationSolve(NetworKit::Graph const &graph) {
        if (graph.numberOfNodes() != graph.upperNodeIdBound()) {
            auto nodeIdMap = NetworKit::GraphTools::getContinuousNodeIds(graph);
            auto compactG = NetworKit::GraphTools::getCompactedGraph(graph, nodeIdMap);
            auto solver = mcp::EnumerationSolverNk(&compactG);
            solver.solve();
            offset_ += solver.bestValue();
        } else {
            auto solver = mcp::EnumerationSolverNk(&graph);
            solver.solve();
            offset_ += solver.bestValue();
        }
        return true;
    }

    double gurobiSolve(NetworKit::Graph const &graph, double timelimit = 3600) {
        if (graph.numberOfNodes() != graph.upperNodeIdBound()) {
            auto nodeIdMap = NetworKit::GraphTools::getContinuousNodeIds(graph);
            auto compactG = NetworKit::GraphTools::getCompactedGraph(graph, nodeIdMap);
            auto solver = GurobiQuadraticSolver(compactG, timelimit, 53115, 4);
            solver.run();
            return solver.getOptimalSolutionValue();
        } else {
            auto solver = GurobiQuadraticSolver(graph, timelimit, 53115, 4);
            solver.run();
            return solver.getOptimalSolutionValue();
        }
    }

    double solutionValue() {
        double solutionValue = offset_;

        for (auto &g : openGraphs_) {
            solutionValue += gurobiSolve(g.graph);
        }

        for (auto &g : finishedGraphs_) {
            solutionValue += gurobiSolve(g);
        }

        return solutionValue;
    }

    NetworKit::Graph const *graph_;
    std::deque<QueueElement> openGraphs_;
    std::vector<NetworKit::Graph> finishedGraphs_;
    bool hasRun_ = false;

    double offset_ = 0.0;

    // config
    int maxConnectivityToUse_;
    int numTrianglesToUse_;
    int cliqueLevel_;
    bool useStableReduce_;
    unsigned int enumThreshold_;
    bool requireFullSolve_;

    // stats
    nlohmann::ordered_json stats_;
    std::chrono::microseconds solvingTime_ = std::chrono::microseconds(0);
};

} // namespace mcp
#endif // MCP_DATA_REDUCTION_HPP
