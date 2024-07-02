#ifndef SMS_NODE_SEPARATOR_HPP
#define SMS_NODE_SEPARATOR_HPP

#include "networkit/graph/Graph.hpp"
#include "networkit/flow/EdmondsKarp.hpp"

class NodeSeparator {

public:
    explicit NodeSeparator(NetworKit::Graph &graph) : graph_(graph), flowLimit_(
            static_cast<NetworKit::edgeweight>(graph_.numberOfNodes()) - 1) {}

    bool protoDelete(NetworKit::node s, unsigned int maxSepSize, unsigned int threshold) {
        initFlowGraph();

        // phase 1
        for (NetworKit::node u: graph_.nodeRange()) {
            if (u != s && !graph_.hasEdge(s, u)) {
                std::vector<int8_t> nodeType;
                auto sepSize = stSep(s, u, nodeType);

                if (sepSize <= maxSepSize) {
                    std::vector<NetworKit::node> leftSide;
                    std::vector<NetworKit::node> rightSide;
                    std::vector<NetworKit::node> separator;
                    for (NetworKit::node curNode: graph_.nodeRange()) {
                        if (nodeType[curNode] == 2) {
                            leftSide.push_back(curNode);
                        } else if (nodeType[curNode] == 1 && curNode != s) {
                            separator.push_back(curNode);
                        } else if (nodeType[curNode] == 0) {
                            rightSide.push_back(curNode);
                        }
                    }

                    assert(separator.size() == sepSize);
                    assert(leftSide.size() > 0);
                    assert(rightSide.size() > 0);

                    if (leftSide.size() < threshold) {
                        std::cout << "Sep size " << sepSize;
                        for (auto u: separator) { std::cout << " ";}

                        std::cout << ". Deleting " << leftSide.size() << " nodes." << std::endl;
                        for (auto leftNode: leftSide) {
                            graph_.removeNode(leftNode);
                        }
                        return true;
                    } else if (rightSide.size() < threshold) {
                        std::cout << "Sep size " << sepSize;
                        std::cout << ". Deleting " << rightSide.size() << " nodes." << std::endl;
                        for (auto rightNode: rightSide) {
                            graph_.removeNode(rightNode);
                        }
                        return true;
                    }
                }
            }
        }

        return false;
    }

    void separator(NetworKit::node s) {
        initFlowGraph();
        std::vector<std::vector<int8_t>> allSeparator;

        // phase 1
        for (NetworKit::node u: graph_.nodeRange()) {
            if (u != s && !graph_.hasEdge(s, u)) {
                std::vector<int8_t> nodeType;
                auto sepSize = stSep(s, u, nodeType);
                if (sepSize < flowLimit_) {
                    allSeparator.emplace_back(std::move(nodeType));
                }
            }
        }

        // phase 2
        for (NetworKit::node u: graph_.neighborRange(s)) {
            for (NetworKit::node v: graph_.neighborRange(s)) {
                if (u < v && !graph_.hasEdge(u, v)) {
                    std::vector<int8_t> nodeType;
                    auto sepSize = stSep(u, v, nodeType);
                    if (sepSize < flowLimit_) {
                        allSeparator.emplace_back(std::move(nodeType));
                    }
                }
            }
        }

        std::cout << "All separators: " << allSeparator.size() << std::endl;

        hasRun_ = true;
    }

    // nodeType is output parameter, sorry
    NetworKit::edgeweight stSep(NetworKit::node s, NetworKit::node t, std::vector<int8_t> &nodeType) {
        // flow from right node to left node
        NetworKit::EdmondsKarp flow(flowGraph_, 2 * s + 1, 2 * t);
        // ToDo: Add an upper bound for the flow (as a parameter to this function) and use it
        flow.run();

        //std::cout << "MaxFlow is: " << flow.getMaxFlow() << std::endl;

        // an edge is part of the cut, if

        if (flow.getMaxFlow() < flowLimit_) {
            nodeType.resize(graph_.upperNodeIdBound(), 0);
            std::fill(nodeType.begin(), nodeType.end(), 0);
            for (auto curNode: flow.getSourceSet()) {
                nodeType[flowToOrig(curNode)]++;
            }
        }

        int leftSideSize = 0;
        int rightSideSize = 0;
        for (NetworKit::node u: graph_.nodeRange()) {
            //std::cout << "(" << i << " : " << int(nodeType[i]) << "), ";
            if (nodeType[u] == 2) {
                leftSideSize++;
            } else if (nodeType[u] == 0) {
                rightSideSize++;
            }
        }
        //std::cout << "left side size is " << leftSideSize << ", right side size is " << rightSideSize << std::endl;

        return flow.getMaxFlow();
    }

    void initFlowGraph() {
        flowGraph_ = NetworKit::Graph(2 * graph_.upperNodeIdBound(), true, true);

        for (NetworKit::node u: graph_.nodeRange()) {
            auto start = 2 * u;
            auto end = start + 1;
            flowGraph_.addEdge(start, end, 1);
        }

        for (auto e: graph_.edgeRange()) {
            auto u = e.u;
            auto v = e.v;
            auto leftU = 2 * u;
            auto rightU = leftU + 1;
            auto leftV = 2 * v;
            auto rightV = leftV + 1;

            flowGraph_.addEdge(rightU, leftV, flowLimit_);
            flowGraph_.addEdge(rightV, leftU, flowLimit_);
        }

        flowGraph_.indexEdges();
    }

    inline static NetworKit::node flowToOrig(NetworKit::node u) {
        return u >> 1;
    }

private:
    bool hasRun_ = false;
    NetworKit::Graph &graph_;
    NetworKit::Graph flowGraph_;
    NetworKit::edgeweight flowLimit_;
};

#endif //SMS_NODE_SEPARATOR_HPP
