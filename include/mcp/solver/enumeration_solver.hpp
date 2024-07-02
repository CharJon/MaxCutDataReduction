#ifndef SMS_ENUMERATION_SOLVER_HPP
#define SMS_ENUMERATION_SOLVER_HPP

#include <vector>

#include "networkit/graph/Graph.hpp"
#include "ogdf/basic/Graph.h"
#include "ogdf/graphalg/steiner_tree/EdgeWeightedGraph.h"

namespace mcp {

/*
 * A class to solve max cut via full enumeration of all solutions.
 * Only works for instances with less than 32 vertices.
 */
class EnumerationSolverNk {
    using partition_t = uint32_t;

public:
    explicit EnumerationSolverNk(const NetworKit::Graph *g) { graph_ = g; }

    /***
     *
     * @param g Graph to solve max cut on
     * @param fixedNodes Vector of (0-indexed) node ids which should be fixed to a certain side
     * @param fixedValues Vector of fixed assignment
     */
    EnumerationSolverNk(const NetworKit::Graph *g, const std::vector<NetworKit::node> &fixedNodes,
                        const std::vector<bool> &fixedValues) {
        assert(fixedNodes.size() == fixedValues.size());

        graph_ = g;

        fixedValues_ = 0;
        for (unsigned int i = 0; i < fixedValues.size(); i++) {
            fixedValues_ |= (fixedValues[i] << fixedNodes[i]);
        }
        fixedMask_ = 0;
        for (auto u : fixedNodes) {
            fixedMask_ |= (1UL << u);
        }
    }

    bool isApplicable() const;

    void solve();

    void solve(const std::vector<NetworKit::node> &fixedNodes, const std::vector<u_int8_t> &fixedValues);

    double bestValue() const;

    std::vector<uint8_t> bestPartitionVec() const;

    NetworKit::edgeweight getCutValue(partition_t) const;

    NetworKit::edgeweight getCutValueLambda(partition_t cut) const;

private:
    NetworKit::Graph const *graph_;
    partition_t bestPartition_ = 0;
    double bestValue_ = 0.;

    partition_t fixedValues_ = 0;
    partition_t fixedMask_ = 0;

    bool hasRun_ = false;

private:
    void solveNoFixed();

    void solveFixed();
};

} // namespace mcp
#endif // SMS_ENUMERATION_SOLVER_HPP
