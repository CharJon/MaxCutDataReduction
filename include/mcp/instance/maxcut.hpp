#ifndef MCP_MAXCUT_HPP
#define MCP_MAXCUT_HPP

#include "networkit/graph/Graph.hpp"

namespace mcp {

/***
 * Represents a MaxCut instance.
 * The graph is compacted, i.e. the node ids are consecutive form [0,...,n-1].
 * The original node ids can be retrieved using getOriginalNode.
 */
class MaxCut {

public:
    explicit MaxCut(const NetworKit::Graph &g, double scalingFactor = 1, double offset = 0);

    unsigned int getNumberOfNodes() const { return graph_.numberOfNodes(); }

    unsigned int getNumberOfEdges() const { return graph_.numberOfEdges(); }

    double getScalingFactor() const { return scalingFactor_; }

    double getOffset() const { return offset_; }

    /***
     * @param solVector Partitioning of nodes.
     * @return Solution value of partition ignoring scaling factor and offset.
     */
    double getSolutionValue(const std::vector<uint8_t> &solVector) const;

    /***
     * @param solVector Partitioning of nodes.
     * @return Solution value of partition including scaling factor and offset.
     */
    double getOriginalSolutionValue(const std::vector<uint8_t> &solVector) const;

    const NetworKit::Graph &getGraph() const { return graph_; }

    NetworKit::node getOriginalNode(NetworKit::node newNode) const { return newToOriginalNode_.at(newNode); }

    NetworKit::node getNewNode(NetworKit::node originalNode) const { return originalToNewNode_.at(originalNode); }

private:
    double scalingFactor_;
    double offset_;
    NetworKit::Graph graph_;
    std::unordered_map<NetworKit::node, NetworKit::node> originalToNewNode_;
    std::unordered_map<NetworKit::node, NetworKit::node> newToOriginalNode_;
};

} // namespace mcp

#endif // MCP_MAXCUT_HPP
