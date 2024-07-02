#ifndef MCP_WEIGHT_STABLE_REDUCER_HPP
#define MCP_WEIGHT_STABLE_REDUCER_HPP

#include <optional>

#include "networkit/clique/MaximalCliques.hpp"
#include "networkit/graph/Graph.hpp"
#include "nlohmann/json.hpp"

#include "mcp/auxiliary/active_elements.hpp"

namespace mcp {

struct RemovableClique {
public:
    std::vector<NetworKit::node> internalNodes;
    std::vector<NetworKit::node> externalNodes;

    size_t size() const { return internalNodes.size() + externalNodes.size(); }
};

using ThreePath = std::array<NetworKit::node, 4>;

/*
 * Data reduction algorithm which does introduce new edge weights for graphs where all edges have the same weight.
 */
class WeightStableReducer {
public:
    /*
     * fullSeparated = Our technique
     * removeSeparated = Any technique of the above
     */
    explicit WeightStableReducer(NetworKit::Graph &g, bool fullSeparatedClique = true,
                                 bool removeSeparatedCliques = true)
        : graph_(g),
          threePathActivity_(graph_.upperNodeIdBound()),
          cliqueActivity_(graph_.upperNodeIdBound()),
          nearCliqueActivity_(graph_.upperNodeIdBound()),
          markings_(g.upperNodeIdBound(), 0),
          fullSeparatedClique_(fullSeparatedClique),
          removeSeparatedCliques_(removeSeparatedCliques) {
        for (auto v : graph_.nodeRange()) {
            threePathActivity_.activate(v);
            cliqueActivity_.activate(v);
            nearCliqueActivity_.activate(v);
        }
    }

    void run();

    int getOffset() const;

    bool hasRun() const;

    bool isClique(const NetworKit::Graph &graph, const RemovableClique &ret) const;

    bool neighborhoodIsClique(NetworKit::node v, bool degreesChecked = false);

    nlohmann::ordered_json getStats() const;

    void sumStats(nlohmann::ordered_json &stats);

    void report(std::ostream &out) const;

    std::optional<RemovableClique> getOneRemovableClique();

    std::optional<RemovableClique> testRemovableNearClique(NetworKit::node u);

private:
    NetworKit::Graph &graph_;
    bool hasRun_ = false;
    int offset_ = 0;

    // active queues
    ActiveElementsStack<NetworKit::node> threePathActivity_;
    ActiveElementsStack<NetworKit::node> cliqueActivity_;
    ActiveElementsStack<NetworKit::node> nearCliqueActivity_;

    // intermediate
    std::vector<uint8_t> markings_;

    // config
    bool fullSeparatedClique_;
    bool removeSeparatedCliques_ = true;

    // stats
    unsigned int numThreePathsContractions_ = 0;
    unsigned int numCliqueRemoval_ = 0;
    unsigned int numCliqueNodesRemoved_ = 0;
    unsigned int numNearCliqueRemoval_ = 0;
    unsigned int numNearCliqueNodesRemoved_ = 0;

    int removeSeparatedCliques();

    void removeClique(const RemovableClique &clique);

    std::optional<ThreePath> testForInducedThreePath();

    void contractThreePath(const ThreePath &indThreePath);

    void markCliqueExternal(NetworKit::node u);

    void activateForAll(NetworKit::node u) {
        threePathActivity_.activate(u);
        cliqueActivity_.activate(u);
        nearCliqueActivity_.activate(u);
    }

    void deactivateForAll(NetworKit::node u) {
        threePathActivity_.deactivate(u);
        cliqueActivity_.deactivate(u);
        nearCliqueActivity_.deactivate(u);
    }

    NetworKit::node partnerNode(NetworKit::node u, NetworKit::node commonNeighbor);
    std::optional<RemovableClique> testRemovableClique(NetworKit::node u);
};

} // namespace mcp

#endif // MCP_WEIGHT_STABLE_REDUCER_HPP
