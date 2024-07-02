#include "mcp/instance/maxcut.hpp"
#include "networkit/graph/GraphTools.hpp"

namespace mcp {

MaxCut::MaxCut(const NetworKit::Graph &g, double scalingFactor, double offset)
    : scalingFactor_(scalingFactor), offset_(offset) {

    NetworKit::node c = 0;
    for (auto v : g.nodeRange()) {
        originalToNewNode_[v] = c;
        newToOriginalNode_[c] = v;
        c++;
    }

    graph_ = NetworKit::GraphTools::getCompactedGraph(g, originalToNewNode_);
    assert(graph_.upperNodeIdBound() == graph_.numberOfNodes());
}

double MaxCut::getSolutionValue(const std::vector<uint8_t> &solVector) const {
    assert(std::all_of(solVector.begin(), solVector.end(), [](int i) { return i == 0 || i == 1; }));
    double res = 0;
    for (auto e : graph_.edgeWeightRange()) {
        res += (solVector[e.u] ^ solVector[e.v]) * e.weight;
    }
    return res;
}

double MaxCut::getOriginalSolutionValue(const std::vector<uint8_t> &solVector) const {
    return getSolutionValue(solVector) * scalingFactor_ + offset_;
}

} // namespace mcp