#include "mcp/solver/burer_solver.hpp"

#include "networkit/auxiliary/Timer.hpp"

mqlib::MaxCutInstance mcp::nkToMqlib(const NetworKit::Graph &g) {
    std::vector<mqlib::Instance::InstanceTuple> mqlibEdgelist;
    for (auto e : g.edgeWeightRange()) {
        mqlibEdgelist.emplace_back(std::make_pair(e.u + 1, e.v + 1), e.weight);
    }
    return mqlib::MaxCutInstance(mqlibEdgelist, static_cast<int>(g.numberOfNodes()));
}

double mcp::BurerMqlibSolver::run() {
    auto t = Aux::Timer();
    t.start();
    run_ = true;
    mqlib::Burer2002 mqHeur(instance_, timelimit_, false, nullptr);
    std::cout << mqHeur.History() << std::endl;
    elapsed_ = t.elapsed();
    bestValue_ = mqHeur.get_best();
    timeToBestSeconds_ = mqHeur.get_best_time();
    pastSolutionValues_ = mqHeur.get_past_solution_values();
    pastSolutionTimes_ = mqHeur.get_past_solution_times();

    return bestValue_;
}

mcp::BurerMqlibSolver::BurerMqlibSolver(const NetworKit::Graph &graph, double timelimit)
    : timelimit_(timelimit), instance_(nkToMqlib(graph)) {}
