#ifndef MCP_SPQR_HPP
#define MCP_SPQR_HPP

#include <optional>

#include "networkit/graph/GraphTools.hpp"
#include "ogdf/basic/Graph.h"
#include "ogdf/decomposition/StaticSPQRTree.h"

namespace mcp {

struct LeafParent {
    ogdf::node leaf;
    ogdf::node parent;
    ogdf::edge edge;
};

class SPQR {
public:
    SPQR(NetworKit::Graph *nkG, unsigned int enumThreshold = 10, bool innerDataReduction = false);

    void run();

    double offset() const { return offset_; }

    unsigned int numLeafContractions() const { return numLeafContractions_; }

    unsigned int numRemovedNodes() const { return numRemovedNodes_; }

    /***
     * Returns the nodes for which neighborhoods changed.
     */
    std::vector<NetworKit::node> activatedNodes() const;

private:
    /*
     * Returns true, if a leaf got processed.
     */
    bool tryToContractLeaf(const ogdf::StaticSPQRTree &spqr, std::vector<uint8_t> &spqrNodeStatus);

    std::optional<LeafParent> contractibleSpqrTreeEdge(const ogdf::StaticSPQRTree &spqr,
                                                       const std::vector<uint8_t> &spqrNodeStatus);

    int numberOfComponents() const { return numberOfComponents_; }

    void spqrTreeReport(const ogdf::StaticSPQRTree &spqr) const;

    NetworKit::Graph *nkG_;
    ogdf::Graph ogdfG_;

    std::unordered_map<NetworKit::node, ogdf::node> nkToOgdfNodeMap_;
    std::unordered_map<int, NetworKit::node> ogdfToNkNodeMap_;
    std::vector<uint8_t > nodeStatus_;

    int numberOfComponents_;
    double offset_ = 0.;

    // config
    unsigned int enumThreshold_ = 2;
    bool innerDataReduction_ = false; // shall data reduction be used for inner nodes of the SPQR tree?

    // stats
    unsigned int numLeafContractions_ = 0;
    unsigned int numRemovedNodes_ = 0;

    void removeTCC(std::vector<NetworKit::node> &nkNodesInSkeleton, NetworKit::node sep1orig, NetworKit::node sep2orig,
                   double bestWhenSame, double bestWhenDiff);
};

} // namespace mcp

#endif // MCP_SPQR_HPP
