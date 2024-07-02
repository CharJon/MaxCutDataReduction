#ifndef MCP_GUROBI_QUADRATIC_SOLVER_HPP
#define MCP_GUROBI_QUADRATIC_SOLVER_HPP

#include "gurobi_c++.h"
#include <string>

#include "networkit/graph/Graph.hpp"
#include "nlohmann/json.hpp"

class GurobiQuadraticSolver {
public:
    GurobiQuadraticSolver(const NetworKit::Graph &graph, double timelimit, int seed = 0, int threads = 1);

    bool hasRun() const;

    double run(const std::string &stats = "", const std::string &logGurobi = "");

    double getRuntimeSeconds() const;

    bool isOptimal() const;

    double getOptimalSolutionValue() const;

    nlohmann::json getStats() const;

private:
    const NetworKit::Graph &graph_;

    bool run_ = false;
    int seed_;
    double timelimit_;
    int threads_;
    GRBEnv env_;
    GRBModel model_;
};

#endif // MCP_GUROBI_QUADRATIC_SOLVER_HPP
