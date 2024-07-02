#ifndef MCP_UNIT_WEIGHT_SOLVER_HPP
#define MCP_UNIT_WEIGHT_SOLVER_HPP

#include <optional>

#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

namespace mcp {

class UnitWeightSolver {
public:
    explicit UnitWeightSolver(const NetworKit::Graph &g,  double unitWeight, bool tryFullKhalfhalfSearch = true);

    bool isApplicable();

    void solve();

    double bestValue() const;

    const std::vector<uint8_t> &bestPartitionVec() const;

private:
    const NetworKit::Graph &graph_;
    bool hasRun_ = false;
    bool applicabilityChecked_ = false;
    bool applicable_ = false;

    // config
    double unitWeight_ = 1;
    bool tryFullKhalfhalfSearch_;

    // data
    std::vector<uint8_t> bestPartition_;
    double bestPartitionValue_ = 0;

    std::optional<std::vector<uint8_t>> getBipartition(const NetworKit::Graph &g) const;
};

} // namespace mcp

#endif // MCP_UNIT_WEIGHT_SOLVER_HPP
