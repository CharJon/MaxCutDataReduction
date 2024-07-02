#include "mcp/prepro/similar_subgraph_reduction.hpp"

namespace mcp {

NetworKit::WeightedEdge edgeToSecondSmallestNeighbor(NetworKit::Graph &g, NetworKit::node u) {
    assert(g.degree(u) > 1);
    NetworKit::node smallestNeighbor = g.getIthNeighbor(u, 0);
    NetworKit::edgeweight smallestNeighborWeight = g.getIthNeighborWeight(u, 0);
    NetworKit::node secondSmallestNeighbor = g.getIthNeighbor(u, 1);
    NetworKit::edgeweight secondSmallestNeighborWeight = g.getIthNeighborWeight(u, 1);
    if (smallestNeighbor > secondSmallestNeighbor) {
        std::swap(smallestNeighbor, secondSmallestNeighbor);
        std::swap(smallestNeighborWeight, secondSmallestNeighborWeight);
    }

    for (NetworKit::count i = 2; i < g.degree(u); i++) {
        auto neighbor = g.getIthNeighbor(u, i);
        auto neighborWeight = g.getIthNeighborWeight(u, i);
        if (neighbor < smallestNeighbor) {
            secondSmallestNeighbor = smallestNeighbor;
            secondSmallestNeighborWeight = smallestNeighborWeight;
            smallestNeighbor = neighbor;
            smallestNeighborWeight = neighborWeight;
        } else if (neighbor < secondSmallestNeighbor) {
            secondSmallestNeighbor = neighbor;
            secondSmallestNeighborWeight = neighborWeight;
        }
    }

    return {u, secondSmallestNeighbor, secondSmallestNeighborWeight};
}

auto SimilarSubgraphReducer::itInBucket(unsigned int bucketNum, unsigned int elementNum) {
    auto it = hashBuckets_.begin(bucketNum);
    // MightDo: is there an easier way, to increase iterator by j?
    for (size_t k = 0; k < elementNum; k++) {
        ++it;
    }
    return it;
}

void SimilarSubgraphReducer::updateEdgeWeight(NetworKit::node currentNode, NetworKit::node keptNeighbor,
                                              NetworKit::node lostNeighbor, NetworKit::edgeweight newWeight) {
    assert(keptNeighbor < lostNeighbor);
    hashBuckets_.erase(currentNode);
    graph_.setWeight(keptNeighbor, currentNode, newWeight);
    nodeInformation_[currentNode].degree -= 1;
    nodeInformation_[currentNode].neighborhoodHash ^= lostNeighbor;
    if (nodeInformation_[currentNode].smallestNeighbor == keptNeighbor)
        nodeInformation_[currentNode].weightSmallestNeighbor = newWeight;
    hashBuckets_.insert(currentNode);
}

void SimilarSubgraphReducer::fillNodeInformation() {
    nodeInformation_.resize(graph_.upperNodeIdBound());
    for (auto v : graph_.nodeRange()) {
        size_t neighborhoodHash = 0;
        if (graph_.degree(v) == 0)
            continue;
        NetworKit::node smallestNeighbor = graph_.getIthNeighbor(v, 0);
        NetworKit::edgeweight smallestNeighborWeight = graph_.getIthNeighborWeight(v, 0);
        for (auto weightedNeighbor : graph_.weightNeighborRange(v)) {
            auto &[neighbor, neighborWeight] = weightedNeighbor;
            neighborhoodHash ^= neighbor;
            if (neighbor < smallestNeighbor) {
                smallestNeighbor = neighbor;
                smallestNeighborWeight = neighborWeight;
            }
        }
        nodeInformation_[v] = {v, graph_.degree(v), neighborhoodHash, smallestNeighbor, smallestNeighborWeight};
    }
}

void SimilarSubgraphReducer::run() {
    while (!activeNodes_.empty()) {
        NetworKit::node u = activeNodes_.popBack();
        if ((graph_.hasNode(u)) && (graph_.degree(u) != 0)) {
            // check bucket of node
            auto nodePair = notConnectedSameNeighborhood(u);
            if (nodePair.has_value()) {
                numSimilarNodesNoEdgeDeleted_ += 1;
                mergeNodesNoEdge(nodePair->at(0), nodePair->at(1));
            } else {
                auto connectedNodePair = connectedSameNeighborhood(u);
                if (connectedNodePair.has_value()) {
                    auto v1 = connectedNodePair.value()[0];
                    auto v2 = connectedNodePair.value()[1];
                    auto type = connectedNodePair.value()[2];
                    switch (type) {
                        case 1:
                        case 2:
                            mergeNodesWithEdge(v1, v2);
                            break;
                        default:
                            throw std::runtime_error("Unknown type.");
                    }
                }
            }
        }

        hasRun_ = true;
    }
}

void SimilarSubgraphReducer::fillHashBuckets() {
    for (auto v : graph_.nodeRange()) {
        hashBuckets_.insert(v);
    }
}

std::optional<std::array<NetworKit::node, 2>> SimilarSubgraphReducer::notConnectedSameNeighborhood() {
    // loop over all buckets unordered_set
    // if (lastBucket >= hashBuckets_.bucket_count())
    unsigned int lastBucket = 0;
    for (; lastBucket < hashBuckets_.bucket_count(); lastBucket++) {
        // list all combinations of elements in this bucket
        auto bucketSize = hashBuckets_.bucket_size(lastBucket);
        if (bucketSize >= 2) {
            for (size_t j = 0; j < bucketSize - 1; j++) {
                auto it = itInBucket(lastBucket, j);
                auto u = *it;
                ++it;
                for (; it != hashBuckets_.end(lastBucket); ++it) {
                    auto v = *it;
                    double alpha = haveSameWeightedNeighborhood(u, v);
                    if (alpha != 0) {
                        return {{u, v}};
                    }
                }
            }
        }
    }
    return {};
}

std::optional<std::array<NetworKit::node, 2>> SimilarSubgraphReducer::notConnectedSameNeighborhood(NetworKit::node u) {
    auto bucketIndex = hashBuckets_.bucket(u);
    // list all combinations of elements in this bucket
    auto bucketSize = hashBuckets_.bucket_size(bucketIndex);
    if (bucketSize >= 2) {
        for (size_t j = 0; j < bucketSize - 1; j++) {
            auto it = itInBucket(bucketIndex, j);
            auto u = *it;
            ++it;
            for (; it != hashBuckets_.end(bucketIndex); ++it) {
                auto v = *it;
                double alpha = haveSameWeightedNeighborhood(u, v);
                if (alpha != 0) {
                    return {{u, v}};
                }
            }
        }
    }
    return {};
}

std::optional<std::array<NetworKit::node, 3>> SimilarSubgraphReducer::connectedSameNeighborhood(NetworKit::node u) {
    for (auto [v, weight] : graph_.weightNeighborRange(u)) {
        if (!activeNodes_.isActive(v)) {
            auto alpha = haveSameWeightedNeighborhood(u, v, true);
            if ((alpha > 0) && (weight < 0)) {
                return {{u, v, 1}};
            } else if ((alpha < 0) && (weight > 0)) {
                return {{u, v, 2}};
            } else if (false && (alpha == 1) && (graph_.degree(u) % 2 == 0)) {
                auto twinWeight = areTwins(u, v);
                if (twinWeight >= weight)
                    return {{u, v, 3}};
            }
            // MightDo: Negative weight twin (alpha == -1)
        }
    }

    return {};
}

// misusing third entry of array to encode type of similar nodes
std::optional<std::array<NetworKit::node, 3>> SimilarSubgraphReducer::connectedSameNeighborhood() {
    for (auto e : graph_.edgeWeightRange()) {
        auto u = e.u;
        auto v = e.v;
        auto weight = e.weight;
        auto alpha = haveSameWeightedNeighborhood(u, v, true);
        if ((alpha > 0) && (weight < 0)) {
            return {{u, v, 1}};
        } else if ((alpha < 0) && (weight > 0)) {
            return {{u, v, 2}};
        } else if (false && (alpha == 1) && (graph_.degree(u) % 2 == 0)) {
            auto twinWeight = areTwins(u, v);
            if (twinWeight >= weight)
                return {{u, v, 3}};
        }
        // MightDo: Negative weight twin (alpha == -1)
    }

    return {};
}

void SimilarSubgraphReducer::mergeNodesNoEdge(NetworKit::node u, NetworKit::node v) {
    assert(!graph_.hasEdge(u, v));

    if (u > v)
        std::swap(u, v); // keep the smaller id, so make u the smaller one
    // hash values get invalidated, so we need to remove and reinsert
    hashBuckets_.erase(u);
    hashBuckets_.erase(v);

    // merge nodes
    auto alpha = nodeInformation_[v].weightSmallestNeighbor / nodeInformation_[u].weightSmallestNeighbor;
    auto alphaAbs = std::abs(alpha);
    double vNeighborSum = 0;
    for (unsigned int i = 0; i < graph_.degree(u); i++) {
        NetworKit::node neighbor = graph_.getIthNeighbor(u, i);
        NetworKit::edgeweight neighborWeight = graph_.getIthNeighborWeight(u, i);
        vNeighborSum += alpha * neighborWeight;
        graph_.setWeight(u, neighbor, neighborWeight * (1 + alphaAbs));

        hashBuckets_.erase(neighbor);
        nodeInformation_[neighbor].degree -= 1;
        nodeInformation_[neighbor].neighborhoodHash ^= v;
        if (nodeInformation_[neighbor].smallestNeighbor == u)
            nodeInformation_[neighbor].weightSmallestNeighbor = neighborWeight * (1 + alphaAbs);
        hashBuckets_.insert(neighbor);
        activeNodes_.activate(neighbor);
    }
    activeNodes_.activate(u);
    assert(alpha != 0);
    double integralPart = 0;
    if (modf(vNeighborSum, &integralPart) != 0) {
        std::cout << "vNeighborSum: " << vNeighborSum << " alpha: " << alpha << " num neighbors: " << graph_.degree(u)
                  << std::endl;
    }
    offset_ += alpha < 0 ? vNeighborSum : 0;
    activeNodes_.deactivate(v);
    graph_.removeNode(v);

    nodeInformation_[u].weightSmallestNeighbor *= (1 + alphaAbs);
    hashBuckets_.insert(u);
}

void SimilarSubgraphReducer::mergeNodesWithEdge(NetworKit::node u, NetworKit::node v) {
    assert(graph_.hasEdge(u, v));
    assert(graph_.degree(u) == graph_.degree(v));

    numSimilarNodesWithEdgeDeleted_ += 1;

    NetworKit::edgeweight connectingEdgeWeight = graph_.weight(u, v);
    hashBuckets_.erase(u);
    hashBuckets_.erase(v);

    if (graph_.degree(u) == 1) {
        offset_ += std::max(0., connectingEdgeWeight);
        graph_.removeNode(u);
        activeNodes_.deactivate(u);
        graph_.removeNode(v);
        activeNodes_.deactivate(v);
        return;
    }

    // delete edge between u and v and update node information
    // update node information of u
    nodeInformation_[u].degree -= 1;
    nodeInformation_[u].neighborhoodHash ^= v;
    if (nodeInformation_[u].smallestNeighbor == v) {
        auto edgeToSecondSmallest = edgeToSecondSmallestNeighbor(graph_, u);
        nodeInformation_[u].smallestNeighbor = edgeToSecondSmallest.v; // v is target node, so correct here!
        nodeInformation_[u].weightSmallestNeighbor = edgeToSecondSmallest.weight;
    }
    // update node information of v
    nodeInformation_[v].degree -= 1;
    nodeInformation_[v].neighborhoodHash ^= u;
    if (nodeInformation_[v].smallestNeighbor == u) {
        auto edgeToSecondSmallest = edgeToSecondSmallestNeighbor(graph_, v);
        nodeInformation_[v].smallestNeighbor = edgeToSecondSmallest.v; // v is target node, so correct here!
        nodeInformation_[v].weightSmallestNeighbor = edgeToSecondSmallest.weight;
    }
    graph_.removeEdge(u, v);

    hashBuckets_.insert(u);
    hashBuckets_.insert(v);

    mergeNodesNoEdge(u, v);

    if (connectingEdgeWeight > 0) {
        assert(connectingEdgeWeight >= 0);
        std::cout << connectingEdgeWeight << std::endl;
        offset_ += connectingEdgeWeight;
    }
}

void SimilarSubgraphReducer::mergeTwins(NetworKit::node u, NetworKit::node v) {
    assert(graph_.hasEdge(u, v));
    assert(graph_.degree(u) % 2 == 0);
    assert(graph_.degree(v) % 2 == 0);

    numTwiNodesDeleted_++;

    if (u > v)
        std::swap(u, v); // keep the smaller id, so make u the smaller one
    // hash values get invalidated, so we need to remove and reinsert
    hashBuckets_.erase(u);
    hashBuckets_.erase(v);
    graph_.removeNode(v);
    activeNodes_.deactivate(v);

    if (graph_.degree(u) == 0)
        return;

    NetworKit::node smallestNeighbor = graph_.getIthNeighbor(u, 0);
    NetworKit::edgeweight smallestNeighborWeight = graph_.getIthNeighborWeight(u, 0);
    for (unsigned int i = 0; i < graph_.degree(u); i++) {
        auto currentNeighbor = graph_.getIthNeighbor(u, i);
        auto neighborWeight = graph_.getIthNeighborWeight(u, i);
        updateEdgeWeight(currentNeighbor, u, v, 2 * neighborWeight);
        if (currentNeighbor < smallestNeighbor) {
            smallestNeighbor = currentNeighbor;
            smallestNeighborWeight = 2 * neighborWeight;
        }
    }

    nodeInformation_[u].degree -= 1;
    nodeInformation_[u].neighborhoodHash ^= v;
    // if (nodeInformation_[u].smallestNeighbor == v) this update is necessary, else it does not hurt
    nodeInformation_[u].smallestNeighbor = smallestNeighbor;
    nodeInformation_[u].weightSmallestNeighbor = smallestNeighborWeight;
    hashBuckets_.insert(u);
}

void SimilarSubgraphReducer::resetMarksOfNeighbors(node u) {
    for (auto neighbor : graph_.neighborRange(u)) {
        marks_[neighbor] = 0;
    }
}

bool SimilarSubgraphReducer::haveSameNeighborhood(NetworKit::node u, NetworKit::node v) {
    assert(std::all_of(marks_.begin(), marks_.end(), [](auto &m) { return m == 0; }));
    assert(graph_.hasNode(u));
    assert(graph_.hasNode(v));
    assert(u != v);

    if (graph_.degree(u) != graph_.degree(v)) {
        return false;
    }
    for (auto neighbor : graph_.neighborRange(u)) {
        marks_[neighbor] = 1;
    }

    bool haveSameNeighborhood = true;
    for (auto neighbor : graph_.neighborRange(v)) {
        if (marks_[neighbor] == 0) {
            haveSameNeighborhood = false;
            break;
        }
    }

    resetMarksOfNeighbors(u);
    return haveSameNeighborhood;
}

double SimilarSubgraphReducer::haveSameWeightedNeighborhood(NetworKit::node u, NetworKit::node v,
                                                            bool ignoreConnection) {
    assert(std::all_of(marks_.begin(), marks_.end(), [](auto &m) { return m == 0; }));
    assert(graph_.hasNode(u));
    assert(graph_.hasNode(v));
    assert(u != v);

    if ((graph_.degree(u) != graph_.degree(v))
        || (ignoreConnection
            && ((nodeInformation_[u].neighborhoodHash ^ v) != (nodeInformation_[v].neighborhoodHash ^ u)))) {
        return 0;
    }

    if (graph_.degree(u) == 1) {
        return 1;
    }

    double alphaU = nodeInformation_[u].weightSmallestNeighbor;
    if (nodeInformation_[u].smallestNeighbor == v) {
        auto edgeToSecondSmallest = edgeToSecondSmallestNeighbor(graph_, u);
        alphaU = edgeToSecondSmallest.weight;
    }

    for (auto weightedNeighbor : graph_.weightNeighborRange(u)) {
        auto &[neighbor, neighborWeight] = weightedNeighbor;
        marks_[neighbor] = neighborWeight / alphaU;
    }
    if (ignoreConnection)
        marks_[v] = 0;

    double alphaV = nodeInformation_[v].weightSmallestNeighbor;
    if (nodeInformation_[v].smallestNeighbor == u) {
        auto edgeToSecondSmallest = edgeToSecondSmallestNeighbor(graph_, v);
        alphaV = edgeToSecondSmallest.weight;
    }

    bool haveSameNeighborhood = true;
    for (auto weightedNeighbor : graph_.weightNeighborRange(v)) {
        auto &[neighbor, neighborWeight] = weightedNeighbor;
        bool ignoreNeighbor = (neighbor == u) && (ignoreConnection);
        if ((not ignoreNeighbor) && (marks_[neighbor] != (neighborWeight / alphaV))) {
            haveSameNeighborhood = false;
            break;
        }
    }

    resetMarksOfNeighbors(u);
    if (haveSameNeighborhood)
        return alphaV / alphaU;
    else
        return 0;
}

double SimilarSubgraphReducer::areTwins(NetworKit::node u, NetworKit::node v) {
    assert(std::all_of(marks_.begin(), marks_.end(), [](auto &m) { return m == 0; }));
    assert(graph_.hasNode(u));
    assert(graph_.hasNode(v));
    assert(u != v);

    if ((graph_.degree(u) != graph_.degree(v)) || (graph_.degree(u) == 0) || (graph_.degree(v) == 0)) {
        return 0.0;
    }

    // check if all neighbors of u have the same weight and safe weight in marks
    auto min = std::numeric_limits<double>::max();
    auto max = std::numeric_limits<double>::min();
    for (auto weightedNeighbor : graph_.weightNeighborRange(u)) {
        auto &[neighbor, neighborWeight] = weightedNeighbor;
        if (neighbor != v) {
            marks_[neighbor] = neighborWeight;
            min = std::min(min, neighborWeight);
            max = std::max(max, neighborWeight);
        }
    }
    if (min != max) {
        resetMarksOfNeighbors(u);
        return 0.0;
    }

    // check if v has the same neighbors as u and if all edges have the same weight
    bool haveSameNeighborhood = true;
    for (auto weightedNeighbor : graph_.weightNeighborRange(v)) {
        auto &[neighbor, neighborWeight] = weightedNeighbor;
        if ((neighbor != u) && (marks_[neighbor] != neighborWeight)) {
            haveSameNeighborhood = false;
            break;
        }
    }

    resetMarksOfNeighbors(u);
    if (haveSameNeighborhood)
        return min; // same as max

    return 0.0;
}

void SimilarSubgraphReducer::removeNode(node u) {
    graph_.removeNode(u);
    activeNodes_.deactivate(u);
}

} // namespace mcp