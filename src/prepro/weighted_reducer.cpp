#include "mcp/prepro/weighted_reducer.hpp"

#include <array>

namespace mcp {

void WeightedReducer::run() {
    hasRun_ = true;
    auto start = std::chrono::high_resolution_clock::now();

    while ((!degreeBuckets_[0].empty()) || (!degreeBuckets_[1].empty()) || (!degreeBuckets_[2].empty())
           || (!degreeBuckets_[3].empty() && useThreeSep_) || (!activeNodes_.empty())) {

        if (!degreeBuckets_[0].empty()) {
            auto v = degreeBuckets_[0].back();
            removeNode(v);
            numDegreeZeroContractions_ += 1;
        } else if (!degreeBuckets_[1].empty()) {
            auto v = degreeBuckets_[1].back();
            reduceDegreeOneNode(v);
        } else if (!degreeBuckets_[2].empty()) {
            auto v = degreeBuckets_[2].back();
            reduceDegreeTwoNode(v);
        } else if (!activeNodes_.empty()) {
            nonDegreeReduction();
        } else if (useThreeSep_ && (!degreeBuckets_[3].empty())) {
            auto v = degreeBuckets_[3].back();
            reduceDegreeThreeNode(v);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    runtime_ = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
}

void WeightedReducer::nonDegreeReduction() {
    auto v = popActiveNode();
    if (tryDominatingEdge(v))
        return;
    if ((trianglesToUse_ != 0) && tryTriangleContraction(v))
        return;
    //  try same neighborhood
}

void WeightedReducer::reduceDegreeOneNode(NetworKit::node v) {
    assert(graph_.degree(v) == 1);
    offset_ += std::max(0., graph_.getIthNeighborWeight(v, 0));
    removeNode(v);
    numDegreeOneContractions_ += 1;
}

void WeightedReducer::reduceDegreeTwoNode(NetworKit::node v) {
    assert(graph_.degree(v) == 2);
    numDegreeTwoContractions_ += 1;

    auto neighbor1 = graph_.getIthNeighbor(v, 0);
    auto neighbor1weight = graph_.getIthNeighborWeight(v, 0);
    auto neighbor2 = graph_.getIthNeighbor(v, 1);
    auto neighbor2weight = graph_.getIthNeighborWeight(v, 1);

    const NetworKit::edgeweight whenInSame = std::max(0., neighbor1weight + neighbor2weight);
    const NetworKit::edgeweight whenInDifferent = std::max(neighbor1weight, neighbor2weight);

    NetworKit::edgeweight newWeight = whenInDifferent - whenInSame;
    if (graph_.hasEdge(neighbor1, neighbor2)) {
        newWeight += graph_.weight(neighbor1, neighbor2);
        removeEdge(neighbor1, neighbor2);
    }

    if (newWeight != 0) {
        addEdge(neighbor1, neighbor2, newWeight);
    }

    offset_ += whenInSame;
    removeNode(v);
}

void WeightedReducer::reduceDegreeThreeNode(NetworKit::node v) {
    const auto neighbor1 = graph_.getIthNeighbor(v, 0);
    const auto weightTo1 = graph_.getIthNeighborWeight(v, 0);
    const auto neighbor2 = graph_.getIthNeighbor(v, 1);
    const auto weightTo2 = graph_.getIthNeighborWeight(v, 1);
    const auto neighbor3 = graph_.getIthNeighbor(v, 2);
    const auto weightTo3 = graph_.getIthNeighborWeight(v, 2);

    NetworKit::edgeweight weight12 = 0;
    if (graph_.hasEdge(neighbor1, neighbor2)) {
        weight12 = graph_.weight(neighbor1, neighbor2);
        removeEdge(neighbor1, neighbor2);
    }

    NetworKit::edgeweight weight13 = 0;
    if (graph_.hasEdge(neighbor1, neighbor3)) {
        weight13 = graph_.weight(neighbor1, neighbor3);
        removeEdge(neighbor1, neighbor3);
    }

    NetworKit::edgeweight weight23 = 0;
    if (graph_.hasEdge(neighbor2, neighbor3)) {
        weight23 = graph_.weight(neighbor2, neighbor3);
        removeEdge(neighbor2, neighbor3);
    }

    const auto optNone = std::max(0., weightTo1 + weightTo2 + weightTo3);
    const auto optNot12 = std::max(weightTo1 + weightTo2, weightTo3) + weight13 + weight23;
    const auto optNot13 = std::max(weightTo2, weightTo1 + weightTo3) + weight12 + weight23;
    const auto optNot23 = std::max(weightTo1, weightTo2 + weightTo3) + weight12 + weight13;

    auto newWeight12 = (optNot13 + optNot23 - optNone - optNot12) / 2;
    if (newWeight12 != 0)
        addEdge(neighbor1, neighbor2, newWeight12);

    auto newWeight13 = (optNot12 + optNot23 - optNone - optNot13) / 2;
    if (newWeight13 != 0)
        addEdge(neighbor1, neighbor3, newWeight13);

    auto newWeight23 = (optNot12 + optNot13 - optNone - optNot23) / 2;
    if (newWeight23 != 0)
        addEdge(neighbor2, neighbor3, newWeight23);

    offset_ += optNone;
    removeNode(v);
    numDegreeThreeContractions_ += 1;
}

void WeightedReducer::contractNodes(NetworKit::node u, NetworKit::node v, bool toSamePartition) {
    assert(graph_.hasNode(u));
    assert(graph_.hasNode(v));

    auto nodeToDelete = graph_.degree(u) < graph_.degree(v) ? u : v; // delete node with smaller degree
    auto nodeToKeep = nodeToDelete == u ? v : u;

    if (!toSamePartition)
        offset_ += graph_.weightedDegree(nodeToDelete);

    for (auto [currentNeighbor, currentWeight] : graph_.weightNeighborRange(nodeToDelete)) {
        assert(graph_.hasNode(currentNeighbor));
        assert(currentNeighbor != nodeToDelete); // no self loop

        if (currentNeighbor != nodeToKeep) {
            NetworKit::edgeweight oldEdgeWeight = 0.;

            if (graph_.hasEdge(currentNeighbor, nodeToKeep)) {
                oldEdgeWeight = graph_.weight(currentNeighbor, nodeToKeep);
                removeEdge(currentNeighbor, nodeToKeep);
            }

            if (toSamePartition && ((oldEdgeWeight + currentWeight) != 0)) {
                addEdge(nodeToKeep, currentNeighbor, (oldEdgeWeight + currentWeight));
            } else if (!toSamePartition && ((oldEdgeWeight - currentWeight) != 0)) {
                addEdge(nodeToKeep, currentNeighbor, (oldEdgeWeight - currentWeight));
            }
        }
    }

    removeNode(nodeToDelete);
}

bool WeightedReducer::tryDominatingEdge(NetworKit::node u) {
    assert(graph_.degree(u) >= 1); // Node has to have at least one neighbor

    NetworKit::edgeweight heaviestWeight = graph_.getIthNeighborWeight(u, 0);
    NetworKit::edgeweight absSum = std::abs(heaviestWeight);
    NetworKit::node neighborWithHighestAbsValue = graph_.getIthNeighbor(u, 0);
    assert(neighborWithHighestAbsValue != NetworKit::none);

    for (unsigned int i = 1; i < graph_.degree(u); i++) {
        const auto currentNeighbor = graph_.getIthNeighbor(u, i);
        assert(currentNeighbor != NetworKit::none);
        const auto currentEdgeWeight = graph_.getIthNeighborWeight(u, i);
        assert(currentEdgeWeight != 0);

        const auto currentAbsEdgeWeight = std::abs(currentEdgeWeight);
        absSum += currentAbsEdgeWeight;
        if (currentAbsEdgeWeight > std::abs(heaviestWeight)) {
            heaviestWeight = currentEdgeWeight;
            neighborWithHighestAbsValue = currentNeighbor;
        }
    }

    if (2 * std::abs(heaviestWeight) >= absSum) {
        contractNodes(u, neighborWithHighestAbsValue, heaviestWeight < 0);
        numHeavyEdgeContractions_ += 1;
        return true;
    }

    return false;
}

bool WeightedReducer::tryTriangleContraction(NetworKit::node u) {
    int maxTriangleCnt = 5;
    int cnt = 0;
    // find triangle
    for (unsigned int i = 0; i < graph_.degree(u); i++) {
        NetworKit::node v = graph_.getIthNeighbor(u, i);
        // We check a triangle only when the last node within it is checked
        if (nodeIsActive(v) == 1)
            continue;
        for (unsigned int j = i + 1; j < graph_.degree(u); j++) {
            if (cnt > maxTriangleCnt)
                return false;
            NetworKit::node w = graph_.getIthNeighbor(u, j);
            if ((!nodeIsActive(w)) && graph_.hasEdge(v, w)) {
                std::array<NetworKit::node, 3> triangle = {u, v, w};
                ++cnt;
                if (tryTriangleContraction(triangle))
                    return true;
            }
        }
    }
    return false;
}

NetworKit::edgeweight absOutWeight(const NetworKit::Graph &graph, NetworKit::node u) {
    NetworKit::edgeweight sum = 0;

    for (auto [neighbor, weight] : graph.weightNeighborRange(u)) {
        sum += std::abs(weight);
    }

    return sum;
}

bool canContractThreeNegative(NetworKit::edgeweight w12, NetworKit::edgeweight w13, NetworKit::edgeweight w23,
                              NetworKit::edgeweight absOutV1, NetworKit::edgeweight absOutV2,
                              NetworKit::edgeweight absOutV3) {
    assert((w12 < 0) && (w13 < 0) && (w23 < 0));

    // e_{1,2} and e_{1,3} are cut
    auto lhs1 = -w12 - w13;
    bool conditionCut1 = lhs1 >= absOutV1 - std::abs(w12) - std::abs(w13);
    bool conditionCut23 = lhs1 >= absOutV2 - std::abs(w12) - std::abs(w23) + absOutV3 - std::abs(w13) - std::abs(w23);
    // e_{1,2} and e_{2,3} are cut
    auto lhs2 = -w12 - w23;
    bool conditionCut2 = lhs2 >= absOutV2 - std::abs(w12) - std::abs(w23);
    bool conditionCut13 = lhs2 >= absOutV1 - std::abs(w12) - std::abs(w13) + absOutV3 - std::abs(w13) - std::abs(w23);

    return (conditionCut1 || conditionCut23) && (conditionCut2 || conditionCut13);
}

// ToDo (performance): Passing weighted triangles would be faster here
bool WeightedReducer::tryTriangleContraction(std::array<NetworKit::node, 3> triangle) {
    auto [v1, v2, v3] = triangle;
    assert(graph_.hasEdge(v1, v2) && graph_.hasEdge(v1, v3) && graph_.hasEdge(v2, v3));

    auto w1 = graph_.weight(v1, v2);
    auto w2 = graph_.weight(v1, v3);
    auto w3 = graph_.weight(v2, v3);

    int numberOfNegativeEdges = (w1 < 0) + (w2 < 0) + (w3 < 0);

    if ((numberOfNegativeEdges == 3) && (trianglesToUse_ & 4u)) {
        NetworKit::edgeweight absOutV1 = absOutWeight(graph_, v1);
        NetworKit::edgeweight absOutV2 = absOutWeight(graph_, v2);
        NetworKit::edgeweight absOutV3 = absOutWeight(graph_, v3);

        if (canContractThreeNegative(w1, w2, w3, absOutV1, absOutV2, absOutV3)) {
            contractNodes(v1, v2, true);
            numTriangleOneContractions_ += 1;
            return true;
        } else if (canContractThreeNegative(w2, w3, w1, absOutV3, absOutV1, absOutV2)) {
            contractNodes(v1, v3, true);
            numTriangleOneContractions_ += 1;
            return true;
        } else if (canContractThreeNegative(w3, w1, w2, absOutV2, absOutV3, absOutV1)) {
            contractNodes(v2, v3, true);
            numTriangleOneContractions_ += 1;
            return true;
        }
    } else if (numberOfNegativeEdges == 1) {
        // make sure that w1 is the one negative edge and w2 and w3 are positive
        if (w2 < 0) {
            std::swap(w1, w2);
            std::swap(v2, v3);
        } else if (w3 < 0) {
            std::swap(w1, w3);
            std::swap(v1, v3);
        }

        NetworKit::edgeweight absOutV1 = absOutWeight(graph_, v1);
        NetworKit::edgeweight absOutV2 = absOutWeight(graph_, v2);
        NetworKit::edgeweight absOutV3 = absOutWeight(graph_, v3);

        if ((trianglesToUse_ & 1u) && (canContractThreeNegative(w1, -w2, -w3, absOutV1, absOutV2, absOutV3))) {
            contractNodes(v1, v2, true);
            numTriangleThreeContractions_ += 1;
            return true;
        } else if ((trianglesToUse_ & 2u) && (canContractThreeNegative(-w2, -w3, w1, absOutV3, absOutV1, absOutV2))) {
            contractNodes(v1, v3, false);
            numTriangleTwoContractions_ += 1;
            return true;
        } else if ((trianglesToUse_ & 2u) && (canContractThreeNegative(-w3, w1, -w2, absOutV2, absOutV3, absOutV1))) {
            contractNodes(v2, v3, false);
            numTriangleTwoContractions_ += 1;
            return true;
        }
    }

    return false;
}

} // namespace mcp