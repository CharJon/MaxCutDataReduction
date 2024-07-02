#ifndef MCP_WEIGHTED_REDUCER_HPP
#define MCP_WEIGHTED_REDUCER_HPP

#include <chrono>

#include "networkit/graph/Graph.hpp"
#include "nlohmann/json.hpp"

namespace mcp {
class WeightedReducer {

public:
    explicit WeightedReducer(NetworKit::Graph &g, bool useThreeSep = true, unsigned int trianglesToUse = 7)
        : graph_(g),
          nodePosInBucket_(g.upperNodeIdBound(), 0),
          nodeIsActive_(g.upperNodeIdBound(), 2),
          nodePosInActiveNodes_(g.upperNodeIdBound(), 0),
          trianglesToUse_(trianglesToUse),
          useThreeSep_(useThreeSep) {
        for (auto u : graph_.nodeRange()) {

            // degree buckets
            auto degree = graph_.degree(u);
            if (degree < 4) {
                nodePosInBucket_[u] = degreeBuckets_[degree].size();
                degreeBuckets_[degree].push_back(u);
            }

            // activity
            markNodeActive(u);
        }
    }

    WeightedReducer(NetworKit::Graph &g, const std::vector<NetworKit::node> &activeNodes, bool useThreeSep = true,
                    unsigned int trianglesToUse = 7)
        : graph_(g),
          nodePosInBucket_(g.upperNodeIdBound(), 0),
          nodeIsActive_(g.upperNodeIdBound(), 2),
          nodePosInActiveNodes_(g.upperNodeIdBound(), 0),
          trianglesToUse_(trianglesToUse),
          useThreeSep_(useThreeSep) {
        for (auto u : graph_.nodeRange()) {

            // degree buckets
            auto degree = graph_.degree(u);
            if (degree < 4) {
                nodePosInBucket_[u] = degreeBuckets_[degree].size();
                degreeBuckets_[degree].push_back(u);
            }
        }

        for (auto u : activeNodes) {
            if (graph_.hasNode(u))
                markNodeActive(u);
        }
    }

    void run();

    bool hasRun() const { return hasRun_; }

    double getOffset() const {
        if (!hasRun_)
            throw std::runtime_error("WeightedReducer has not run yet.");
        return offset_;
    }

    int numTotalReductions() const {
        return numDegreeZeroContractions_ + numDegreeOneContractions_ + numDegreeTwoContractions_
               + numDegreeThreeContractions_ + numHeavyEdgeContractions_ + numTriangleOneContractions_
               + numTriangleTwoContractions_ + numTriangleThreeContractions_;
    }

    auto getTotalRuntime() const { return runtime_; }

    double getTotalRuntimeSeconds() const {
        std::chrono::duration<double> totalRuntimeSeconds = std::chrono::duration_cast<std::chrono::seconds>(runtime_);
        return totalRuntimeSeconds.count();
    }

    nlohmann::ordered_json getStats() const {
        nlohmann::ordered_json stats;
        stats["numDegreeZeroContractions"] = numDegreeZeroContractions_;
        stats["numDegreeOneContractions"] = numDegreeOneContractions_;
        stats["numDegreeTwoContractions"] = numDegreeTwoContractions_;
        stats["numDegreeThreeContractions"] = numDegreeThreeContractions_;
        stats["numHeavyEdgeContractions"] = numHeavyEdgeContractions_;
        stats["numTriangleOneContractions"] = numTriangleOneContractions_;
        stats["numTriangleTwoContractions"] = numTriangleTwoContractions_;
        stats["numTriangleThreeContractions"] = numTriangleThreeContractions_;
        stats["numTotalReductions"] = numTotalReductions();
        stats["totalRuntimeSeconds"] = getTotalRuntimeSeconds();
        return stats;
    }

    void sumStats(nlohmann::ordered_json &stats) const {
        stats["numDegreeZeroContractions"] =
            static_cast<unsigned int>(stats["numDegreeZeroContractions"]) + numDegreeZeroContractions_;
        stats["numDegreeOneContractions"] =
            static_cast<unsigned int>(stats["numDegreeOneContractions"]) + numDegreeOneContractions_;
        stats["numDegreeTwoContractions"] =
            static_cast<unsigned int>(stats["numDegreeTwoContractions"]) + numDegreeTwoContractions_;
        stats["numDegreeThreeContractions"] =
            static_cast<unsigned int>(stats["numDegreeThreeContractions"]) + numDegreeThreeContractions_;
        stats["numHeavyEdgeContractions"] =
            static_cast<unsigned int>(stats["numHeavyEdgeContractions"]) + numHeavyEdgeContractions_;
        stats["numTriangleOneContractions"] =
            static_cast<unsigned int>(stats["numTriangleOneContractions"]) + numTriangleOneContractions_;
        stats["numTriangleTwoContractions"] =
            static_cast<unsigned int>(stats["numTriangleTwoContractions"]) + numTriangleTwoContractions_;
        stats["numTriangleThreeContractions"] =
            static_cast<unsigned int>(stats["numTriangleThreeContractions"]) + numTriangleThreeContractions_;
        stats["totalRuntimeSeconds"] = static_cast<double>(stats["totalRuntimeSeconds"]) + getTotalRuntimeSeconds();
    }

    void report(std::ostream &out) const {
        out << "WeightedReducer stats:" << std::endl;
        out << getStats().dump(4) << std::endl;
    }

private:
    NetworKit::Graph &graph_;
    double offset_ = 0;

    bool hasRun_ = false;

    // data structures
    std::array<std::vector<NetworKit::node>, 4> degreeBuckets_;
    std::vector<size_t> nodePosInBucket_;

    std::vector<uint8_t> nodeIsActive_;
    std::vector<NetworKit::node> activeNodes_;
    std::vector<size_t> nodePosInActiveNodes_;

    // config
    unsigned int trianglesToUse_ = 7; // binary flag, 7 is all, 4 is rule one only etc.
    bool useThreeSep_;

    // stats
    std::chrono::microseconds runtime_ = std::chrono::microseconds(0);
    int numDegreeZeroContractions_ = 0;
    int numDegreeOneContractions_ = 0;
    int numDegreeTwoContractions_ = 0;
    int numDegreeThreeContractions_ = 0;
    int numHeavyEdgeContractions_ = 0;
    int numTriangleOneContractions_ = 0;
    int numTriangleTwoContractions_ = 0;
    int numTriangleThreeContractions_ = 0;

    bool nodeIsActive(NetworKit::node u) { return nodeIsActive_[u] == 1; }

    void markNodeActive(NetworKit::node u) {
        if (!nodeIsActive(u)) {
            nodeIsActive_[u] = 1;
            nodePosInActiveNodes_[u] = activeNodes_.size();
            activeNodes_.push_back(u);
        }
    }

    NetworKit::node popActiveNode() {
        assert(!activeNodes_.empty());
        auto u = activeNodes_.back();
        assert(nodeIsActive(u));
        assert(activeNodes_[nodePosInActiveNodes_[u]] == u);
        activeNodes_.pop_back();
        nodeIsActive_[u] = 0;
        return u;
    }

    void removeEdge(NetworKit::node u, NetworKit::node v) {
        decreaseDegreeByOne(u);
        decreaseDegreeByOne(v);
        markNodeActive(u);
        markNodeActive(v);
        graph_.removeEdge(u, v);
    }

    void addEdge(NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        assert(weight != 0);
        increaseDegreeByOne(u);
        increaseDegreeByOne(v);
        markNodeActive(u);
        markNodeActive(v);
        graph_.addEdge(u, v, weight);
    }

    void removeNode(NetworKit::node u) {
        removeNodeFromBuckets(u);
        for (auto neighbor : graph_.neighborRange(u)) {
            decreaseDegreeByOne(neighbor);
            markNodeActive(neighbor);
        }

        if (nodeIsActive(u)) {
            deactivateNode(u);
        }

        graph_.removeNode(u);
    }

    void deactivateNode(NetworKit::node u) {
        auto nodePos = nodePosInActiveNodes_[u];
        std::swap(activeNodes_[nodePos], activeNodes_.back());
        nodePosInActiveNodes_[activeNodes_[nodePos]] = nodePos;
        activeNodes_.pop_back();
    }

    void decreaseDegreeByOne(NetworKit::node u) {
        auto degree = graph_.degree(u);
        if (degree == 4) {
            // new bucket
            auto &newBucket = degreeBuckets_[degree - 1];
            nodePosInBucket_[u] = newBucket.size();
            newBucket.push_back(u);
        } else if (degree <= 3) {
            auto &oldBucket = degreeBuckets_[degree];
            auto posInOldBucket = nodePosInBucket_[u];
            assert(oldBucket[posInOldBucket] == u);
            std::swap(oldBucket[posInOldBucket], oldBucket.back());
            nodePosInBucket_[oldBucket[posInOldBucket]] = posInOldBucket;
            oldBucket.pop_back();
            auto &newBucket = degreeBuckets_[degree - 1];
            nodePosInBucket_[u] = newBucket.size();
            newBucket.push_back(u);
        }
    }

    void increaseDegreeByOne(NetworKit::node u) {
        auto degree = graph_.degree(u);
        if (degree <= 3) {
            // old bucket
            auto &oldBucket = degreeBuckets_[degree];
            auto posInOldBucket = nodePosInBucket_[u];
            assert(oldBucket[posInOldBucket] == u);
            std::swap(oldBucket[posInOldBucket], oldBucket.back());
            nodePosInBucket_[oldBucket[posInOldBucket]] = posInOldBucket;
            oldBucket.pop_back();
        }
        if (degree <= 2) {
            // new bucket
            auto &newBucket = degreeBuckets_[degree + 1];
            nodePosInBucket_[u] = newBucket.size();
            newBucket.push_back(u);
        }
    }

    void removeNodeFromBuckets(NetworKit::node u) {
        auto degree = graph_.degree(u);
        if (degree <= 3) {
            // old bucket
            auto &oldBucket = degreeBuckets_[degree];
            auto posInOldBucket = nodePosInBucket_[u];
            assert(oldBucket[posInOldBucket] == u);
            std::swap(oldBucket[posInOldBucket], oldBucket.back());
            nodePosInBucket_[oldBucket[posInOldBucket]] = posInOldBucket;
            oldBucket.pop_back();
        }
    }

    void reduceDegreeOneNode(NetworKit::node v);

    void reduceDegreeTwoNode(NetworKit::node v);

    void reduceDegreeThreeNode(NetworKit::node v);

    void contractNodes(NetworKit::node u, NetworKit::node v, bool toSamePartition);

    void nonDegreeReduction();

    bool tryDominatingEdge(NetworKit::node u);

    bool tryTriangleContraction(NetworKit::node u);

    bool tryTriangleContraction(std::array<NetworKit::node, 3> triangle);
};
} // namespace mcp

#endif // MCP_WEIGHTED_REDUCER_HPP
