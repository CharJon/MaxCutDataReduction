#ifndef MCP_BASE_HPP
#define MCP_BASE_HPP

#include "networkit/graph/Graph.hpp"

namespace mcp {

void removeZeroWeightAndZeroDegree(NetworKit::Graph &g) {
    std::vector<NetworKit::Edge> edgesToRemove;
    for (const auto &e : g.edgeWeightRange()) {
        if (e.weight == 0) {
            edgesToRemove.emplace_back(e.u, e.v);
        }
    }
    for (auto &e : edgesToRemove) {
        g.removeEdge(e.u, e.v);
    }

    std::vector<NetworKit::node> nodesToRemove;
    for (const auto n : g.nodeRange()) {
        if (g.degree(n) == 0) {
            nodesToRemove.push_back(n);
        }
    }
    for (const auto n : nodesToRemove) {
        g.removeNode(n);
    }
}

} // namespace mcp

#endif // MCP_BASE_HPP
