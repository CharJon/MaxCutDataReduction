#include "mcp/solver/enumeration_solver.hpp"

#include <vector>

#include "networkit/graph/Graph.hpp"

namespace mcp {

template <typename T>
inline bool inDifferentPartitions(T cut, unsigned int pos1, unsigned int pos2) {
    return ((cut >> pos1) ^ (cut >> pos2)) & (1UL);
}

bool EnumerationSolverNk::isApplicable() const {
    return (graph_->numberOfNodes() == graph_->upperNodeIdBound()) && (graph_->numberOfNodes() < 32);
}

void EnumerationSolverNk::solve() {
    if (!isApplicable())
        throw std::runtime_error("EnumerationSolver is not applicable to this graph");

    hasRun_ = true;
    if (graph_->numberOfNodes() == 0) {
        bestValue_ = 0;
        bestPartition_ = 0;
    } else if (fixedMask_ == 0) {
        solveNoFixed();
    } else {
        solveFixed();
    }
}

void EnumerationSolverNk::solve(const std::vector<NetworKit::node> &fixedNodes,
                                const std::vector<uint8_t> &fixedValues) {
    assert(fixedNodes.size() == fixedValues.size());
    if (!isApplicable())
        throw std::runtime_error("EnumerationSolver is not applicable to this graph");

    hasRun_ = true;
    fixedValues_ = 0;
    for (unsigned int i = 0; i < fixedValues.size(); i++) {
        fixedValues_ |= (fixedValues[i] << fixedNodes[i]);
    }
    fixedMask_ = 0;
    for (auto u : fixedNodes) {
        fixedMask_ |= (1UL << u);
    }
    solveFixed();
}

void EnumerationSolverNk::solveNoFixed() {
    assert(fixedMask_ == 0);

    partition_t currentPartition = 0UL;

    partition_t best = currentPartition;
    double currentBestValue = getCutValue(currentPartition);

    partition_t upperBoundSolutionVector = 1UL << (graph_->numberOfNodes() - 1);
    while (currentPartition < upperBoundSolutionVector) {
        double currentValue = getCutValueLambda(currentPartition);

        if (currentValue > currentBestValue) {
            best = currentPartition;
            currentBestValue = currentValue;
        }

        currentPartition += 1;
    }

    bestPartition_ = best;
    bestValue_ = currentBestValue;
}

void EnumerationSolverNk::solveFixed() {
    assert(fixedMask_ != 0);

    partition_t currentPartition = fixedValues_;

    partition_t best = currentPartition;
    double currentBestValue = getCutValue(currentPartition);

    currentPartition |= fixedMask_; // set all fixed to one
    currentPartition += 1;

    partition_t upperBoundSolutionVector = 1UL << graph_->numberOfNodes();
    while (currentPartition < upperBoundSolutionVector) {
        currentPartition &= (~fixedMask_);  // reset fixed
        currentPartition |= (fixedValues_); // set fixed to correct values

        double currentValue = getCutValueLambda(currentPartition);

        if (currentValue > currentBestValue) {
            best = currentPartition;
            currentBestValue = currentValue;
        }

        currentPartition |= fixedMask_; // set all fixed to one
        currentPartition += 1;
    }

    bestPartition_ = best;
    bestValue_ = currentBestValue;
}

double EnumerationSolverNk::getCutValue(partition_t cut) const {
    double cutValue = 0;

    for (auto e : graph_->edgeWeightRange()) {
        cutValue += inDifferentPartitions(cut, e.u, e.v) ? e.weight : 0;
    }

    return cutValue;
}

inline double EnumerationSolverNk::getCutValueLambda(partition_t cut) const {
    double cutValue = 0;

    graph_->forEdges([&cutValue, &cut](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        cutValue += inDifferentPartitions(cut, u, v) * weight;
    });

    return cutValue;
}

std::vector<uint8_t> EnumerationSolverNk::bestPartitionVec() const {
    assert(hasRun_);
    auto sol = std::vector<uint8_t>(graph_->numberOfNodes(), false);
    for (auto u : graph_->nodeRange()) {
        sol[u] = (bestPartition_ >> u) & 1UL;
    }
    return sol;
}

double EnumerationSolverNk::bestValue() const {
    assert(hasRun_);
    return bestValue_;
}

} // namespace mcp
