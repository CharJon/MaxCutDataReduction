#include "mcp/prepro/spqr.hpp"

#include "networkit/auxiliary/Timer.hpp"

#include "mcp/solver/enumeration_solver.hpp"

namespace mcp {

SPQR::SPQR(NetworKit::Graph *nkG, unsigned int enumThreshold, bool innerDataReduction)
    : nkG_(nkG),
      nodeStatus_(nkG->upperNodeIdBound(), 0),
      enumThreshold_(enumThreshold),
      innerDataReduction_(innerDataReduction) {
    for (auto u : nkG_->nodeRange()) {
        auto v = ogdfG_.newNode();
        nkToOgdfNodeMap_[u] = v;
        ogdfToNkNodeMap_[v->index()] = u;
    }

    for (auto e : nkG_->edgeRange()) {
        auto u = e.u;
        auto v = e.v;
        assert(ogdfG_.searchEdge(nkToOgdfNodeMap_[u], nkToOgdfNodeMap_[v]) == nullptr);
        ogdfG_.newEdge(nkToOgdfNodeMap_[u], nkToOgdfNodeMap_[v]);
    }
}

void SPQR::run() {
    // start chrono timer here
    auto t = Aux::Timer();
    t.start();

    ogdf::StaticSPQRTree spqr(ogdfG_);
    std::cout << "SPQR tree creation " << t.elapsedTag() << std::endl;
    spqrTreeReport(spqr);

    auto spqrNodeStatus = std::vector<uint8_t>(0);
    for (auto u : spqr.tree().nodes) {
        if (static_cast<unsigned long>(u->index()) >= spqrNodeStatus.size()) {
            spqrNodeStatus.resize(u->index() + 1, 1);
        }
    }
    // assert all element in spqrNodeStatus are 1, 0 means deleted, 2 means can not be solved
    assert(std::all_of(spqrNodeStatus.begin(), spqrNodeStatus.end(), [](auto &a) { return a == 1; }));

    // run leaf contraction
    assert(spqr.tree().numberOfNodes() > 0);
    bool foundOne = spqr.tree().numberOfNodes() > 0;
    int numRounds = 0;
    while (foundOne) {
        numRounds += 1;
        if (numRounds % 100 == 0) {
            std::cout << "Round " << numRounds << std::endl;
        }
        foundOne = tryToContractLeaf(spqr, spqrNodeStatus);
    }

    std::cout << "SPQR reduction " << t.elapsedTag() << std::endl;
}

bool SPQR::tryToContractLeaf(const ogdf::StaticSPQRTree &spqr, std::vector<uint8_t> &spqrNodeStatus) {
    auto leafAndParent = contractibleSpqrTreeEdge(spqr, spqrNodeStatus);
    if (leafAndParent.has_value()) {
        auto leaf = leafAndParent.value().leaf;
        // auto parent = leafParent.value().parent;
        auto &leafSkeleton = spqr.skeleton(leaf);

        if (spqr.typeOf(leafAndParent.value().leaf) == ogdf::SPQRTree::NodeType::PNode) {
            spqrNodeStatus[leaf->index()] = 0;
            return true;
        }

        if (static_cast<unsigned int>(spqr.skeleton(leaf).getGraph().numberOfNodes()) > enumThreshold_) {
            spqrNodeStatus[leaf->index()] = 2;
            return true;
        }

        ogdf::edge edgeInLeafSkeleton;
        if (leaf == leafAndParent.value().edge->source()) {
            edgeInLeafSkeleton = spqr.skeletonEdgeSrc(leafAndParent.value().edge);
        } else {
            edgeInLeafSkeleton = spqr.skeletonEdgeTgt(leafAndParent.value().edge);
        }

        // collect all nodes in skeleton with original id
        std::vector<NetworKit::node> nkNodesInSkeleton;
        auto &currentSkeleton = spqr.skeleton(leaf);
        for (auto v : currentSkeleton.getGraph().nodes) {
            nkNodesInSkeleton.push_back(ogdfToNkNodeMap_[currentSkeleton.original(v)->index()]);
        }

        // extract subgraph based on nodes from original graph
        // ToDo (performance): Theoretically unnecessary copy of graph
        auto subGraph =
            NetworKit::GraphTools::subgraphFromNodes(*nkG_, nkNodesInSkeleton.begin(), nkNodesInSkeleton.end(), false);
        auto continuesNodeIds = NetworKit::GraphTools::getContinuousNodeIds(subGraph);
        auto subGraphCompact = NetworKit::GraphTools::getCompactedGraph(subGraph, continuesNodeIds);

        auto uInLeafSkeleton = leafSkeleton.original(edgeInLeafSkeleton->source());
        auto vInLeafSkeleton = leafSkeleton.original(edgeInLeafSkeleton->target());
        NetworKit::node sep1orig = ogdfToNkNodeMap_[uInLeafSkeleton->index()];
        NetworKit::node sep2orig = ogdfToNkNodeMap_[vInLeafSkeleton->index()];

        if (subGraphCompact.hasEdge(continuesNodeIds.at(sep1orig), continuesNodeIds.at(sep2orig))) {
            subGraphCompact.removeEdge(continuesNodeIds.at(sep1orig), continuesNodeIds.at(sep2orig));
        }

        // if problem is small enough, solve it via enumeration
        if (subGraphCompact.numberOfNodes() <= enumThreshold_) {
            // solve subgraph with two partitions
            auto solver = mcp::EnumerationSolverNk(&subGraphCompact);
            solver.solve({continuesNodeIds.at(sep1orig), continuesNodeIds.at(sep2orig)}, {0, 0});
            auto bestWhenSame = solver.bestValue();
            solver.solve({continuesNodeIds.at(sep1orig), continuesNodeIds.at(sep2orig)}, {0, 1});
            auto bestWhenDiff = solver.bestValue();

            removeTCC(nkNodesInSkeleton, sep1orig, sep2orig, bestWhenSame, bestWhenDiff);
            spqrNodeStatus[leaf->index()] = 0;
        } else {
            spqrNodeStatus[leaf->index()] = 2;
        }
        return true;
    }

    return false;
}

void SPQR::removeTCC(std::vector<NetworKit::node> &nkNodesInSkeleton, NetworKit::node sep1orig,
                     NetworKit::node sep2orig, double bestWhenSame, double bestWhenDiff) {
    numLeafContractions_ += 1;
    numRemovedNodes_ += nkNodesInSkeleton.size() - 2;

    nodeStatus_[sep1orig] = 2;
    nodeStatus_[sep2orig] = 2;

    for (auto nodeToDelete : nkNodesInSkeleton) {
        if (nodeToDelete != sep1orig && nodeToDelete != sep2orig) {
            nodeStatus_[nodeToDelete] = 1;
            nkG_->removeNode(nodeToDelete);
        }
    }

    offset_ += bestWhenSame;
    if (bestWhenDiff - bestWhenSame != 0) {
        if (nkG_->hasEdge(sep1orig, sep2orig)) {
            auto oldWeight = nkG_->weight(sep1orig, sep2orig);
            double newWeight = oldWeight + bestWhenDiff - bestWhenSame;
            nkG_->removeEdge(sep1orig, sep2orig);
            if (newWeight != 0) {
                nkG_->addEdge(sep1orig, sep2orig, newWeight);
            }
        } else {
            nkG_->addEdge(sep1orig, sep2orig, bestWhenDiff - bestWhenSame);
        }
    }
}

std::optional<LeafParent> SPQR::contractibleSpqrTreeEdge(const ogdf::StaticSPQRTree &spqr,
                                                         const std::vector<uint8_t> &spqrNodeStatus) {
    // ToDo (performance): This is a linear search. We could use a degree one bucket.
    for (auto u : spqr.tree().nodes) {
        if (spqrNodeStatus[u->index()] != 1) {
            continue;
        }

        int numRemainingNeighbors = 0;
        ogdf::node parent = nullptr;
        ogdf::edge spqrEdge = nullptr;
        for (auto adj : u->adjEntries) {
            auto v = adj->theEdge()->opposite(u);
            if (spqrNodeStatus[v->index()] > 0) {
                parent = v;
                spqrEdge = adj->theEdge();
                numRemainingNeighbors++;
            }
        }

        if (numRemainingNeighbors == 1) {
            LeafParent leafParent(u, parent, spqrEdge);
            return {leafParent};
        }
    }
    return {};
}

void SPQR::spqrTreeReport(const ogdf::StaticSPQRTree &spqr) const {
    std::cout << "SPQR tree has " << spqr.tree().numberOfNodes() << " nodes "
              << " and " << spqr.tree().numberOfEdges() << " edges." << std::endl;

    for (auto u : spqr.tree().nodes) {
        auto &g = spqr.skeleton(u).getGraph();
        if (static_cast<unsigned int>(g.numberOfNodes()) > enumThreshold_)
            std::cout << "Node " << u << " is of size " << g.numberOfNodes() << "." << std::endl;
    }
}

std::vector<NetworKit::node> SPQR::activatedNodes() const {
    std::vector<NetworKit::node> activatedNodes;
    for (auto u : nkG_->nodeRange()) {
        if (nodeStatus_[u] == 2)
            activatedNodes.push_back(u);
    }

    return activatedNodes;
}

} // namespace mcp