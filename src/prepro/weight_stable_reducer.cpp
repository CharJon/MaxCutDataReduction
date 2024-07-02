#include "mcp/prepro/weight_stable_reducer.hpp"

#include <concepts>
#include <map>
#include <optional>

namespace mcp {

// Divide an integer by two and ceil the result
template <std::integral T>
inline T divTwoCeil(T a) {
    return (a + T{1}) / 2;
}

int cutSizePerfectSplit(unsigned int numNodes) {
    int numSmallerSide = static_cast<int>(numNodes / 2);             // div 2 floor
    int numLargerSide = static_cast<int>(numNodes) - numSmallerSide; // equal or larger by one
    return numSmallerSide * numLargerSide;
}

bool WeightStableReducer::neighborhoodIsClique(NetworKit::node v, bool degreesChecked) {
    assert(graph_.hasNode(v));
    assert(graph_.hasNode(v));
    assert(std::all_of(markings_.begin(), markings_.end(), [](auto i) { return i == 0; }));

    if (not degreesChecked) {
        for (NetworKit::node w : graph_.neighborRange(v)) {
            if (graph_.degree(w) < graph_.degree(v))
                return false;
        }
    }

    bool isGood = true;
    unsigned int cliqueSize = graph_.degree(v) + 1;

    markings_[v] = 1;
    for (NetworKit::node w : graph_.neighborRange(v))
        markings_[w] = 1;

    for (NetworKit::node w : graph_.neighborRange(v)) {
        unsigned int i = 0;
        // Count the number of clique neighbors of w
        for (NetworKit::node x : graph_.neighborRange(w)) {
            i += markings_[x];
        }
        // Check if w neighbors every node in the supposed clique
        if (i != cliqueSize - 1) {
            isGood = false;
            break;
        }
    }

    // Restore the markings
    markings_[v] = 0;
    for (NetworKit::node w : graph_.neighborRange(v))
        markings_[w] = 0;

    return isGood;
}

bool isCliqueC(const NetworKit::Graph &graph, const std::vector<NetworKit::node> &nodes) {
    // faster with markings
    for (unsigned int i = 0; i < nodes.size(); i++) {
        for (unsigned int j = i + 1; j < nodes.size(); j++) {
            if (!graph.hasEdge(nodes[i], nodes[j]))
                return false;
        }
    }
    return true;
}

bool WeightStableReducer::isClique(const NetworKit::Graph &graph, const RemovableClique &ret) const {
    for (NetworKit::node v : ret.externalNodes) {
        for (NetworKit::node w : ret.externalNodes) {
            if (v != w && !graph.hasEdge(v, w))
                return false;
        }
    }
    for (NetworKit::node v : ret.internalNodes) {
        for (NetworKit::node w : ret.internalNodes) {
            if (v != w && !graph.hasEdge(v, w))
                return false;
        }
    }
    for (NetworKit::node v : ret.internalNodes) {
        for (NetworKit::node w : ret.externalNodes) {
            if (v != w && !graph.hasEdge(v, w))
                return false;
        }
    }
    return true;
}

std::optional<RemovableClique> WeightStableReducer::testRemovableClique(NetworKit::node v) {
    auto potentialCliqueSize = graph_.degree(v) + 1;

    bool defNotIt = false;
    // Step one for checking if this is the internal vertex of an acceptable clique:
    // There may not be a neighbor with smaller degree
    unsigned int numLargeNeighbors = 0;
    for (NetworKit::node w : graph_.neighborRange(v)) {
        if (graph_.degree(w) > graph_.degree(v)) {
            numLargeNeighbors++;
            cliqueActivity_.deactivate(w);
        } else if (graph_.degree(w) == graph_.degree(v))
            cliqueActivity_.deactivate(w);
        else {
            // If v has a neighbor with lower degree its definitely not an internal vertex of a clique
            defNotIt = true;
            break;
        }
    }

    // This checks if number of largeNeighbors (which would all be external) > ceil(clique/2)
    if (numLargeNeighbors > divTwoCeil(potentialCliqueSize)) {
        defNotIt = true;
    }

    if (defNotIt) {
        // If v isn't internal to a removable clique then none of its neighbors with degree >= are either
        for (NetworKit::node w : graph_.neighborRange(v)) {
            if (graph_.degree(w) >= graph_.degree(v))
                cliqueActivity_.deactivate(w);
        }
        return {};
    }

    // Now we know for certain that if the neighborhood of v is a clique that it is a removable one
    if (neighborhoodIsClique(v, true)) {
        RemovableClique ret;
        ret.externalNodes.reserve(numLargeNeighbors);
        ret.internalNodes.reserve(graph_.degree(v) + 1 - numLargeNeighbors);
        ret.internalNodes.emplace_back(v);
        for (NetworKit::node w : graph_.neighborRange(v)) {
            if (graph_.degree(w) == graph_.degree(v))
                ret.internalNodes.emplace_back(w);
            else
                ret.externalNodes.emplace_back(w);
        }
        assert(ret.externalNodes.size() <= divTwoCeil(ret.externalNodes.size() + ret.internalNodes.size()));

        return ret;
    }

    // If the neighborhood of v isn't a removable clique (and all vertices have degree >= deg(v)), they are
    // all not internal to a removable clique
    for (NetworKit::node w : graph_.neighborRange(v))
        cliqueActivity_.deactivate(w);
    return {};
}

std::optional<RemovableClique> WeightStableReducer::getOneRemovableClique() {
    while (!cliqueActivity_.empty()) {
        auto v = cliqueActivity_.popBack();
        if (!graph_.hasNode(v)) {
            deactivateForAll(v);
            continue;
        }

        auto res = testRemovableClique(v);
        if (res.has_value()) {
            return res;
        }
    }
    return {};
}

// find a node v that has the same neighborhood as u and share the common neighbor
NetworKit::node WeightStableReducer::partnerNode(NetworKit::node u, NetworKit::node commonNeighbor) {
    assert(std::all_of(markings_.begin(), markings_.end(), [](auto u) { return u == 0; }));
    assert(graph_.hasEdge(u, commonNeighbor));
    for (NetworKit::node v : graph_.neighborRange(commonNeighbor)) {
        if ((v != u) && (graph_.degree(u) == graph_.degree(v))) {
            bool found = true;

            for (auto w : graph_.neighborRange(u)) {
                markings_[w] = 1;
            }

            for (auto w : graph_.neighborRange(v)) {
                found &= (markings_[w] == 1);
            }

            for (auto w : graph_.neighborRange(u)) {
                markings_[w] = 0;
            }

            if (found) {
                return v;
            }
        }
    }

    return NetworKit::none;
}

std::optional<RemovableClique> WeightStableReducer::testRemovableNearClique(NetworKit::node u) {
    if (graph_.degree(u) == 0) {
        deactivateForAll(u);
        graph_.removeNode(u);
        return {};
    }

    NetworKit::node minNeigh = graph_.getIthNeighbor(u, 0);
    NetworKit::count minDeg = graph_.degree(minNeigh);
    NetworKit::count numLargeNeighbors = 0;
    auto potentialCliqueSize = graph_.degree(u) + 2;

    // If any neighbor of v has a lower or equal degree
    // we can immediately say that v is not internal for a near clique
    for (NetworKit::node w : graph_.neighborRange(u)) {
        if (graph_.degree(w) <= graph_.degree(u))
            return {};
        if (graph_.degree(w) >= graph_.degree(u)) {
            if (graph_.degree(w) >= potentialCliqueSize)
                numLargeNeighbors++;
        }
        minDeg = std::min(minDeg, graph_.degree(w));
        if (graph_.degree(w) == minDeg)
            minNeigh = w;
    }

    // Since we are exclusively looking for near-cliques with odd size or at least 3 internal vertices and at most
    // half external vertices we immediately check for this criteria
    if ((potentialCliqueSize % 2 != 0) || potentialCliqueSize - numLargeNeighbors >= 2) {
        // we are allowed to add the edge, this does not guarantee for the clique to be removable
        // here we need to make sure not to many nodes are external
        if (numLargeNeighbors > divTwoCeil(potentialCliqueSize)) {
            return {};
        }

        if (neighborhoodIsClique(u, true)) {
            auto sameNeighborhoodNode = partnerNode(u, minNeigh);
            if (sameNeighborhoodNode != NetworKit::none) {
                assert(!graph_.hasEdge(u, sameNeighborhoodNode));
                graph_.addEdge(u, sameNeighborhoodNode);
                // create and return clique
                std::vector<NetworKit::node> internalNodes = {u};
                std::vector<NetworKit::node> externalNodes = {};
                for (auto v : graph_.neighborRange(u)) {
                    if (graph_.degree(v) > graph_.degree(u)) {
                        externalNodes.push_back(v);
                    } else {
                        internalNodes.push_back(v);
                    }
                }
                return {RemovableClique(std::move(internalNodes), std::move(externalNodes))};
            }
        }

        return {};
    }

    return {};
}

void WeightStableReducer::run() {
    bool didDelete = true;

    while (didDelete) {
        didDelete = false;

        auto indThreePath = testForInducedThreePath();
        if (indThreePath.has_value()) {
            didDelete = true;
            contractThreePath(indThreePath.value());
            numThreePathsContractions_++;
        }

        if (!didDelete) {
            auto clique = getOneRemovableClique();
            if (clique.has_value()) {
                didDelete = true;
                numCliqueRemoval_++;
                numCliqueNodesRemoved_ += clique.value().internalNodes.size();
                removeClique(clique.value());
            }
        }

        if (!didDelete) {
            while (!nearCliqueActivity_.empty()) {
                auto v = nearCliqueActivity_.popBack();
                if (!graph_.hasNode(v))
                    continue;
                auto nearClique = testRemovableNearClique(v);
                if (nearClique.has_value()) {
                    didDelete = true;
                    numNearCliqueRemoval_++;
                    numNearCliqueNodesRemoved_ += nearClique.value().size();
                    removeClique(nearClique.value());
                    break;
                }
            }
        }
    }
    if (removeSeparatedCliques_)
        removeSeparatedCliques();

    hasRun_ = true;
}

void WeightStableReducer::contractThreePath(const ThreePath &indThreePath) {
    assert(graph_.hasEdge(indThreePath[0], indThreePath[1]));
    assert(graph_.hasEdge(indThreePath[1], indThreePath[2]));
    assert(graph_.hasEdge(indThreePath[2], indThreePath[3]));
    assert(!graph_.hasEdge(indThreePath[0], indThreePath[3]));
    graph_.addEdge(indThreePath[0], indThreePath[3], 1.0, false);

    // activate for cliques
    activateForAll(indThreePath[0]);
    activateForAll(indThreePath[3]);
    for (auto u : graph_.neighborRange(indThreePath[0])) {
        activateForAll(u);
    }
    for (auto u : graph_.neighborRange(indThreePath[3])) {
        activateForAll(u);
    }

    graph_.removeNode(indThreePath[1]);
    deactivateForAll(indThreePath[1]);
    graph_.removeNode(indThreePath[2]);
    deactivateForAll(indThreePath[2]);
    offset_ += 2;
}

void WeightStableReducer::removeClique(const RemovableClique &clique) {
    assert(isClique(graph_, clique));
    unsigned int cliqueSize = clique.internalNodes.size() + clique.externalNodes.size();
    unsigned int numCutEdges = cutSizePerfectSplit(cliqueSize);
    offset_ += static_cast<int>(numCutEdges);

    for (auto u : clique.internalNodes) {
        graph_.removeNode(u);
        deactivateForAll(u);
    }
    for (unsigned int i = 0; i < clique.externalNodes.size(); i++) {
        for (unsigned int j = i + 1; j < clique.externalNodes.size(); j++) {
            auto u = clique.externalNodes[i];
            auto v = clique.externalNodes[j];
            assert(graph_.hasEdge(u, v));
            graph_.removeEdge(u, v);
        }
    }
    for (auto u : clique.externalNodes) {
        deactivateForAll(u);

        if (not graph_.hasNode(u)) {
            // skip if u was removed in the meantime
        } else if (graph_.degree(u) == 1) {
            auto currentNode = u;
            while (graph_.degree(currentNode) == 1) {
                auto neighbor = graph_.getIthNeighbor(currentNode, 0);
                graph_.removeNode(currentNode);
                currentNode = neighbor;
                offset_ += 1;
            }
            if (graph_.degree(currentNode) == 0) {
                graph_.removeNode(currentNode);
            } else {
                markCliqueExternal(currentNode);
            }
        } else {
            markCliqueExternal(u);
        }
    }
}

void WeightStableReducer::markCliqueExternal(NetworKit::node u) {
    if (graph_.degree(u) == 2) {
        threePathActivity_.activate(u);
    }
    bool maybeCliqueInternal = true;

    // neighbors of degree 2 get marked as active
    // as they may now be part of an induced three paths
    for (auto neighbor : graph_.neighborRange(u)) {
        if (graph_.degree(neighbor) == 2) {
            threePathActivity_.activate(neighbor);
        }
        // Neighbors with a lower degree are marked as candidates for cliques
        if (graph_.degree(neighbor) < graph_.degree(u)) {
            cliqueActivity_.activate(neighbor);
            nearCliqueActivity_.activate(neighbor);
            maybeCliqueInternal = false;
        }
        // If no neighbor has a lower degree we mark u
        // If u and a neighbor have the same degree they are either both internal to the same clique or neither is,
        // so marking u is enough
        if (maybeCliqueInternal) {
            cliqueActivity_.activate(u);
            nearCliqueActivity_.activate(u);
        }
    }
}

int WeightStableReducer::getOffset() const {
    assert(hasRun_);
    return offset_;
}

bool WeightStableReducer::hasRun() const {
    return hasRun_;
}

std::optional<ThreePath> WeightStableReducer::testForInducedThreePath() {
    while (!threePathActivity_.empty()) {
        auto nextNode = threePathActivity_.popBack();
        if (graph_.hasNode(nextNode) && (graph_.degree(nextNode) == 2)) {
            auto leftNeighbor = graph_.getIthNeighbor(nextNode, 0);
            auto rightNeighbor = graph_.getIthNeighbor(nextNode, 1);

            if (graph_.degree(leftNeighbor) == 2) {
                auto leftLeft = graph_.getIthNeighbor(leftNeighbor, 0);
                auto leftRight = graph_.getIthNeighbor(leftNeighbor, 1);
                auto end = leftLeft == nextNode ? leftRight : leftLeft;
                if ((end != rightNeighbor) && !graph_.hasEdge(end, rightNeighbor)) {
                    return {{rightNeighbor, nextNode, leftNeighbor, end}};
                }
            }

            if (graph_.degree(rightNeighbor) == 2) {
                auto rightLeft = graph_.getIthNeighbor(rightNeighbor, 0);
                auto rightRight = graph_.getIthNeighbor(rightNeighbor, 1);
                auto end = rightLeft == nextNode ? rightRight : rightLeft;
                if ((end != leftNeighbor) && !graph_.hasEdge(end, leftNeighbor)) {
                    return {{leftNeighbor, nextNode, rightNeighbor, end}};
                }
            }
        }
    }
    return {};
}

void WeightStableReducer::report(std::ostream &out) const {
    if (!hasRun_) {
        throw std::runtime_error("WeightStableReducer::report(): Reduction has not been run yet.");
    }

    out << "WeightStableReducer: " << std::endl;
    out << getStats().dump(4) << std::endl;
}

nlohmann::ordered_json WeightStableReducer::getStats() const {
    nlohmann::ordered_json stats;
    stats["numThreePathsContractions"] = numThreePathsContractions_;
    stats["numCliqueRemoval"] = numCliqueRemoval_;
    stats["numCliqueNodesRemoved"] = numCliqueNodesRemoved_;
    stats["numNearCliqueRemoval"] = numNearCliqueRemoval_;
    stats["numNearCliqueNodesRemoved"] = numNearCliqueNodesRemoved_;
    return stats;
}

void WeightStableReducer::sumStats(nlohmann::ordered_json &stats) {
    stats["numThreePathsContractions"] =
        static_cast<unsigned int>(stats["numThreePathsContractions"]) + numThreePathsContractions_;
    stats["numCliqueRemoval"] = static_cast<unsigned int>(stats["numCliqueRemoval"]) + numCliqueRemoval_;
    stats["numCliqueNodesRemoved"] = static_cast<unsigned int>(stats["numCliqueNodesRemoved"]) + numCliqueNodesRemoved_;
    stats["numNearCliqueRemoval"] = static_cast<unsigned int>(stats["numNearCliqueRemoval"]) + numNearCliqueRemoval_;
    stats["numNearCliqueNodesRemoved"] =
        static_cast<unsigned int>(stats["numNearCliqueNodesRemoved"]) + numNearCliqueNodesRemoved_;
}

int WeightStableReducer::removeSeparatedCliques() {
    // hash nodes based on their neighborhood + themselves
    auto buckets = std::map<std::vector<NetworKit::node>, std::vector<NetworKit::node>>();
    // collect buckets as values of hashmap
    for (auto u : graph_.nodeRange()) {
        std::vector<NetworKit::node> neighbors(graph_.neighborRange(u).begin(), graph_.neighborRange(u).end());
        neighbors.push_back(u);
        std::sort(neighbors.begin(), neighbors.end());

        if (buckets.find(neighbors) == buckets.end()) {
            buckets[neighbors] = std::vector<NetworKit::node>{u};
        } else {
            buckets[neighbors].push_back(u);
        }
    }

    int c = 0;
    // loop over each bucket
    // for each one it is guaranteed that the nodes are in the same clique and have the same neighbors
    for (const auto &bucket : buckets) {
        // check total size of neighborhood, which for any node u is deg(u) +1 - bucketSize
        // if neighborhood size is less or equal to bucket size, then the clique is removable
        if ((bucket.second.size() > 1) && (2 * bucket.second.size() + 1 >= bucket.first.size())) {
            assert(isCliqueC(graph_, bucket.second));
            c++;
            if (fullSeparatedClique_) {
                // collect neighborhood of clique
                // mark nodes which are part of the clique
                for (auto u : bucket.second) {
                    markings_[u] = 1;
                }
                // unmarked nodes form neighborhood
                std::vector<NetworKit::node> neighborhood;
                for (auto u : bucket.first) {
                    if (markings_[u] == 0) {
                        neighborhood.push_back(u);
                    }
                }
                // reset marks
                for (auto u : bucket.second) {
                    markings_[u] = 0;
                }

                for (unsigned int i = 0; i < neighborhood.size(); i++) {
                    for (unsigned int j = i + 1; j < neighborhood.size(); j++) {
                        if (graph_.hasEdge(neighborhood[i], neighborhood[j])) {
                            graph_.removeEdge(neighborhood[i], neighborhood[j]);
                        } else {
                            graph_.addEdge(neighborhood[i], neighborhood[j], -1.);
                        }
                    }
                }

                for (auto u : bucket.second) {
                    graph_.removeNode(u);
                }
                offset_ += cutSizePerfectSplit(bucket.first.size());
            } else {
                std::vector<NetworKit::node> remainingCliqueNodes = bucket.second;
                unsigned int curNeighborhoodSize = bucket.first.size() - remainingCliqueNodes.size();

                while (remainingCliqueNodes.size() > curNeighborhoodSize && remainingCliqueNodes.size() > 1) {
                    // remove two nodes from clique
                    auto x1 = remainingCliqueNodes.back();
                    remainingCliqueNodes.pop_back();
                    auto x2 = remainingCliqueNodes.back();
                    remainingCliqueNodes.pop_back();

                    offset_ += static_cast<int>(graph_.degree(x1));
                    graph_.removeNode(x1);
                    graph_.removeNode(x2);
                }

                if (remainingCliqueNodes.size() == curNeighborhoodSize && !remainingCliqueNodes.empty()) {
                    offset_ += static_cast<int>(remainingCliqueNodes.size());
                    auto x = remainingCliqueNodes.back();
                    remainingCliqueNodes.pop_back();
                    graph_.removeNode(x);
                }
            }
        }
    }

    return c;
}

} // namespace mcp
