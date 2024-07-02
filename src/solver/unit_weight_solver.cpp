#include "mcp/solver/unit_weight_solver.hpp"

#include "mcp/auxiliary/graphs.hpp"

namespace mcp {

UnitWeightSolver::UnitWeightSolver(const NetworKit::Graph &g, double unitWeight, bool tryFullKhalfhalfSearch)
    : graph_(g),
      unitWeight_(unitWeight),
      tryFullKhalfhalfSearch_(tryFullKhalfhalfSearch),
      bestPartition_(g.upperNodeIdBound(), 3) {}

bool UnitWeightSolver::isApplicable() {
    applicabilityChecked_ = true;

    if (graph_.numberOfNodes() != graph_.upperNodeIdBound())
        throw std::runtime_error("Not implemented for non compact graphs.");
#ifndef NDEBUG
    for (auto e : graph_.edgeWeightRange()) {
        assert(e.weight == unitWeight_);
    }
#endif

    NetworKit::count n = graph_.numberOfNodes();
    NetworKit::count m = graph_.numberOfEdges();
    NetworKit::count numPossibleEdges = (n * (n - 1)) / 2;
    NetworKit::count numEdgesKayHalfHalf = (n / 2) * ((n + 1) / 2);

    // Complete graphs are easy
    if (m == numPossibleEdges) {
        for (unsigned int i = 0; i < n / 2; i++) {
            bestPartition_[i] = 0;
        }
        for (unsigned int i = n / 2; i < n; i++) {
            bestPartition_[i] = 1;
        }
        bestPartitionValue_ = unitWeight_ * static_cast<double>(numEdgesKayHalfHalf);
        applicable_ = true;
        return true;
    }

    // Bipartite graphs are easy
    auto bipartition = getBipartition(graph_);
    if (bipartition.has_value()) {
        bestPartition_ = bipartition.value();
        bestPartitionValue_ = unitWeight_ * static_cast<double>(graph_.numberOfEdges());
        applicable_ = true;
        return true;
    }

    if ((m >= (numPossibleEdges - n / 2)) || (tryFullKhalfhalfSearch_ && (m >= numEdgesKayHalfHalf))) {
        // Graph may contain a K_{n/2, n/2}, so we get the complement graph
        NetworKit::Graph complement = mcp::unweightedComplementGraph(graph_); //  getComplementGraph(graph_);

        NetworKit::ConnectedComponents components(complement);
        components.run();
        auto connComps = components.getComponents();

        std::vector<unsigned int> componentSizes(components.numberOfComponents());

        NetworKit::count compNum = 0;

        for (const std::vector<NetworKit::node> &comp : connComps) {
            if (comp.size() > (n + 1) / 2) {
                applicable_ = false;
                return false;
            }
            componentSizes[compNum] = comp.size();
            compNum++;
        }

        std::vector<std::vector<bool>> doable(components.numberOfComponents());
        for (auto &i : doable)
            i.resize((n + 1) / 2 + 1, false);

        // We see if the CCs in the complement graph can be divided into equal size groups
        doable[0][0] = true;
        doable[0][componentSizes[0]] = true;

        for (unsigned int i = 1; i < doable.size(); i++) {
            for (unsigned int j = 0; j < (n + 1) / 2 + 1; j++) {
                if (static_cast<int>(j) - static_cast<int>(componentSizes[i]) >= 0)
                    doable[i][j] = doable[i - 1][j] || doable[i - 1][j - componentSizes[i]];
                else
                    doable[i][j] = doable[i - 1][j];
            }
        }

        // Now we can check if the graph contains a K_{n/2, n/2}
        if (doable[componentSizes.size() - 1][(n + 1) / 2]) {
            bestPartitionValue_ = unitWeight_ * static_cast<double>(numEdgesKayHalfHalf);

            NetworKit::count pos = (n + 1) / 2;
            for (NetworKit::count j = componentSizes.size() - 1; j > 0; j--) {
                assert(doable[j][pos]);
                if (doable[j][pos] == doable[j - 1][pos]) {
                    for (NetworKit::node v : connComps[j])
                        bestPartition_[v] = 1;
                } else {
                    for (NetworKit::node v : connComps[j])
                        bestPartition_[v] = 0;
                    pos -= connComps[j].size();
                }
            }

            if (pos == 0)
                for (NetworKit::node v : connComps[0])
                    bestPartition_[v] = 1;
            else
                for (NetworKit::node v : connComps[0])
                    bestPartition_[v] = 0;

            applicable_ = true;
            return true;
        }
    }

    applicable_ = false;
    return false;
}

void UnitWeightSolver::solve() {
    if (not applicabilityChecked_)
        isApplicable();

    if (not applicable_) {
        throw std::runtime_error("Solver is not applicable!");
    }

    // the actual work is done in the applicability check
    hasRun_ = true;
}

double UnitWeightSolver::bestValue() const {
    if (!hasRun_)
        throw std::runtime_error("Call run() first!");
    return bestPartitionValue_;
}

const std::vector<uint8_t> &UnitWeightSolver::bestPartitionVec() const {
    if (!hasRun_)
        throw std::runtime_error("Call run() first!");
    return bestPartition_;
}

std::optional<std::vector<uint8_t>> UnitWeightSolver::getBipartition(const NetworKit::Graph &g) const {
    std::vector<uint8_t> partition(g.upperNodeIdBound(), 3);

    for (NetworKit::node v : g.nodeRange()) {
        if (partition[v] != 3)
            continue;

        partition[v] = 0;

        std::deque queue(1, v);

        while (not queue.empty()) {
            NetworKit::node w = queue.front();
            queue.pop_front();

            assert(partition[w] != 3);

            for (NetworKit::node x : g.neighborRange(w)) {
                if (partition[w] == partition[x]) {
                    return {};
                }

                if (partition[x] == 3) {
                    queue.emplace_back(x);
                    partition[x] = 1 - partition[w];
                }
            }
        }
    }

    return partition;
}

} // namespace mcp
