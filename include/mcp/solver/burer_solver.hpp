#ifndef MCP_BURER_SOLVER_HPP
#define MCP_BURER_SOLVER_HPP

#include "chrono"

#include "mqlib/heuristics/maxcut/burer2002.h"
#include "networkit/graph/Graph.hpp"

namespace mcp {

class BurerMqlibSolver {
public:
    BurerMqlibSolver(const NetworKit::Graph &graph, double timelimit);

    bool hasRun() const { return run_; }

    double run();

    std::chrono::duration<uint64_t, std::milli> elapsed() const {
        if (!run_) {
            throw std::runtime_error("Solver has not been run yet.");
        }
        return elapsed_;
    }

    double bestValue() const {
        if (!run_) {
            throw std::runtime_error("Solver has not been run yet.");
        }
        return bestValue_;
    }

    double timeToBestSeconds() const {
        if (!run_) {
            throw std::runtime_error("Solver has not been run yet.");
        }
        return timeToBestSeconds_;
    }

    const std::vector<double> &pastSolutionValues() const {
        if (!run_) {
            throw std::runtime_error("Solver has not been run yet.");
        }
        return pastSolutionValues_;
    }

    const std::vector<double> &pastSolutionTimes() const {
        if (!run_) {
            throw std::runtime_error("Solver has not been run yet.");
        }
        return pastSolutionTimes_;
    }

private:
    bool run_ = false;
    std::chrono::duration<uint64_t, std::milli> elapsed_;
    double timelimit_;
    mqlib::MaxCutInstance instance_;

    double bestValue_ = 0;
    double timeToBestSeconds_ = 0;
    std::vector<double> pastSolutionValues_;
    std::vector<double> pastSolutionTimes_;
};

mqlib::MaxCutInstance nkToMqlib(const NetworKit::Graph &g);

} // namespace mcp

#endif // MCP_BURER_SOLVER_HPP