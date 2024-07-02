#include "mcp/solver/gurobi_quadratic_solver.hpp"

#include "gurobi_c++.h"
#include <fstream>

#include "nlohmann/json.hpp"

GurobiQuadraticSolver::GurobiQuadraticSolver(const NetworKit::Graph &graph, double timelimit, int seed, int threads)
    : graph_(graph), seed_(seed), timelimit_(timelimit), threads_(threads), env_(), model_(env_) {
    assert(graph_.numberOfNodes() == graph_.upperNodeIdBound());
    assert(graph_.checkConsistency());
}

bool GurobiQuadraticSolver::hasRun() const {
    return run_;
}

double GurobiQuadraticSolver::getRuntimeSeconds() const {
    if (!run_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }
    return model_.get(GRB_DoubleAttr_Runtime);
}

bool GurobiQuadraticSolver::isOptimal() const {
    if (!run_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }
    return model_.get(GRB_IntAttr_Status) == GRB_OPTIMAL;
}

double GurobiQuadraticSolver::getOptimalSolutionValue() const {
    if (!run_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }
    if (model_.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        throw std::runtime_error("Not optimal.");
    }
    return model_.get(GRB_DoubleAttr_ObjVal);
}

nlohmann::json GurobiQuadraticSolver::getStats() const {
    if (!run_) {
        throw std::runtime_error("GurobiQuadraticSolver has not run yet.");
    }
    nlohmann::json j;

    j["cpu_solving_time"] = model_.get(GRB_DoubleAttr_Runtime);
    j["best_solution_value"] = model_.get(GRB_DoubleAttr_ObjVal);
    j["optimal_solution_found"] = model_.get(GRB_IntAttr_Status) == GRB_OPTIMAL;
    j["timelimit_reached"] = model_.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT;
    j["bab_nodes"] = model_.get(GRB_DoubleAttr_NodeCount);
    j["obj_bound_value"] = model_.get(GRB_DoubleAttr_ObjBound);
    j["seed"] = seed_;
    j["timelimit"] = timelimit_;
    j["threads"] = threads_;

    return j;
}

double GurobiQuadraticSolver::run(const std::string &stats, const std::string &logGurobi) {
    run_ = true;
    int numberOfNodes = static_cast<int>(graph_.numberOfNodes());

    try {
        GRBVar *vars = model_.addVars(numberOfNodes, GRB_BINARY);

        if (!logGurobi.empty()) {
            model_.set(GRB_StringParam_LogFile, logGurobi);
            model_.set(GRB_IntParam_LogToConsole, 0);
        }

        model_.set(GRB_IntParam_Seed, seed_);
        model_.set(GRB_DoubleParam_TimeLimit, timelimit_);
        model_.set(GRB_IntParam_Threads, threads_);
        model_.set(GRB_DoubleParam_MIPGap, 1e-6);
        model_.set(GRB_IntParam_Presolve, 2);
        // model_.set(GRB_IntParam_Disconnected, 2);

        GRBQuadExpr quadObjective(0.0);

        for (int i = 0; i < numberOfNodes; i++) {
            quadObjective.addTerm(graph_.weightedDegree(i), vars[i], vars[i]);
            //vars[i].set(GRB_IntAttr_BranchPriority, static_cast<int>(graph_.degree(i)));
        }

        for (int i = 0; i < numberOfNodes; i++) {
            for (int j = i + 1; j < numberOfNodes; j++) {
                auto matrixEntry = 2 * graph_.weight(i, j);
                if (matrixEntry != 0.0) {
                    quadObjective.addTerm(-matrixEntry, vars[i], vars[j]);
                }
            }
        }

        model_.setObjective(quadObjective, GRB_MAXIMIZE);

        model_.optimize();

        if (!stats.empty()) {
            std::ofstream o(stats);
            o << std::setw(4) << getStats() << std::endl;
        }

        run_ = true;

        auto finalSolValue = model_.get(GRB_DoubleAttr_ObjVal);
        delete vars;

        return finalSolValue;

    } catch (GRBException &e) {
        std::cout << "Error code = " << e.getErrorCode() << '\n' << e.getMessage() << std::endl;
        throw e;
    }
}
