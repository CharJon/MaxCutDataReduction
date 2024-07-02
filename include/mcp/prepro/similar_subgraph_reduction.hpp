#ifndef MCP_SIMILAR_SUBGRAPH_REDUCTION_HPP
#define MCP_SIMILAR_SUBGRAPH_REDUCTION_HPP

#include <optional>

#include "networkit/graph/Graph.hpp"

#include "mcp/auxiliary/active_elements.hpp"
#include "mcp/auxiliary/graphs.hpp"

namespace mcp {

class SimilarSubgraphReducer {
public:
    explicit SimilarSubgraphReducer(NetworKit::Graph &g)
        : graph_(g),
          marks_(g.upperNodeIdBound(), 0),
          hasher_(this),
          hashBuckets_(graph_.numberOfNodes(), hasher_),
          activeNodes_(graph_.upperNodeIdBound()) {
        assert(std::none_of(graph_.nodeRange().begin(), graph_.nodeRange().end(),
                            [&g](auto v) { return g.degree(v) == 0; }));
        fillNodeInformation();
        fillHashBuckets();
        for (auto v : graph_.nodeRange()) {
            activeNodes_.activate(v);
        }
    }

    SimilarSubgraphReducer(const SimilarSubgraphReducer &) = delete;

    SimilarSubgraphReducer(SimilarSubgraphReducer &&) = delete;

    void run();

    bool hasRun() const { return hasRun_; }

    double getOffset() const { return offset_; }

    unsigned int numReductions() const {
        return numSimilarNodesNoEdgeDeleted_ + numSimilarNodesWithEdgeDeleted_ + numTwiNodesDeleted_;
    }

    bool haveSameNeighborhood(NetworKit::node u, NetworKit::node v);

    /***
     * Edge between u and v is ignored
     * @return alpha if u and v have the same weighted neighborhood, 0 otherwise
     */
    double haveSameWeightedNeighborhood(NetworKit::node u, NetworKit::node v, bool ignoreConnection = false);

    /***
     * @return 0 if u and v are not twins, the weight of all edges leaving u and v otherwise
     */
    double areTwins(node u, node v);

    void mergeTwins(NetworKit::node u, NetworKit::node v);

    std::optional<std::array<NetworKit::node, 2>> notConnectedSameNeighborhood();

    std::optional<std::array<NetworKit::node, 2>> notConnectedSameNeighborhood(node u);

    std::optional<std::array<NetworKit::node, 3>> connectedSameNeighborhood();

    std::optional<std::array<NetworKit::node, 3>> connectedSameNeighborhood(node u);

    void report(std::ostream &out) const {
        out << "SimilarSubgraphReducer:" << std::endl;
        out << "  numSimilarNodesNoEdgeDeleted: " << numSimilarNodesNoEdgeDeleted_ << std::endl;
        out << "  numSimilarNodesWithEdgeDeleted: " << numSimilarNodesWithEdgeDeleted_ << std::endl;
        out << "  numTwinNodesDeleted: " << numTwiNodesDeleted_ << std::endl;
        out << "  offset: " << offset_ << std::endl;
        // out << "  alphaCases: " << alphaCases_ << std::endl;
    }

private:
    constexpr static uint8_t kINACTIVE = 0;
    constexpr static uint8_t kACTIVE = 1;

    struct NodeInformation {
        NetworKit::node node;
        NetworKit::count degree;
        size_t neighborhoodHash;
        NetworKit::node smallestNeighbor;
        NetworKit::edgeweight weightSmallestNeighbor;
    };

    struct NodeInformationHash {
        explicit NodeInformationHash(SimilarSubgraphReducer *outer) : outer_(outer){};
        size_t operator()(NetworKit::node v) const { return outer_->nodeInformation_[v].neighborhoodHash; }

    private:
        SimilarSubgraphReducer *outer_;
    };

    NetworKit::Graph &graph_;
    bool hasRun_ = false;
    // intermediate data structures
    std::vector<NetworKit::edgeweight> marks_;
    std::vector<NodeInformation> nodeInformation_;
    NodeInformationHash hasher_;
    std::unordered_set<NetworKit::node, NodeInformationHash> hashBuckets_;
    ActiveElementsStack<NetworKit::node> activeNodes_;

    // stats
    double offset_ = 0;
    unsigned int numSimilarNodesNoEdgeDeleted_ = 0;
    unsigned int numSimilarNodesWithEdgeDeleted_ = 0;
    unsigned int numTwiNodesDeleted_ = 0;

    void resetMarksOfNeighbors(node u);
    void fillNodeInformation();
    void fillHashBuckets();
    void mergeNodesNoEdge(node u, node v);
    void mergeNodesWithEdge(node u, node v);

    auto itInBucket(unsigned int bucketNum, unsigned int elementNum);
    void updateEdgeWeight(node currentNode, node keptNeighbor, node lostNeighbor, edgeweight newWeight);

    void removeNode(node u);
};

} // namespace mcp

#endif // MCP_SIMILAR_SUBGRAPH_REDUCTION_HPP
