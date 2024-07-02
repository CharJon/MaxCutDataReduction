#ifndef SMS_BLOCK_CUT_TREE_HPP
#define SMS_BLOCK_CUT_TREE_HPP

#include <map>
#include <unordered_map>
#include <vector>

#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

class BlockCutTree {
public:
    explicit BlockCutTree(const NetworKit::Graph &graph) : graph_(&graph), biconnected_components_(*graph_) {}

    void run() {
        biconnected_components_.run();

        articulation_points_of_component_.resize(biconnected_components_.numberOfComponents());

        for (NetworKit::node u : graph_->nodeRange()) {
            auto &componentsOfU = biconnected_components_.getComponentsOfNode(u);
            if (componentsOfU.size() > 1) {
                articulation_points_.push_back(u);

                for (auto &curComponent : componentsOfU) {
                    articulation_points_of_component_[curComponent].push_back(u);
                }
            }
        }

        components_ = biconnected_components_.getComponents();

        has_run_ = true;
    }

    NetworKit::count numberOfComponents() const { return components_.size(); }

    NetworKit::count numberOfArticulationPoints() const { return articulation_points_.size(); }

    const std::vector<std::vector<NetworKit::node>> &getComponents() const {
        assert(has_run_);
        return components_;
    };

    const std::vector<NetworKit::node> &getArticulationPoints() const {
        assert(has_run_);
        return articulation_points_;
    }

    const std::unordered_set<NetworKit::node> &getComponentsOfNode(NetworKit::node u) const {
        assert(has_run_);
        return biconnected_components_.getComponentsOfNode(u);
    }

    const std::vector<NetworKit::node> &getNodesOfComponent(NetworKit::count component) const {
        assert(has_run_);
        return components_[component];
    }

    const std::vector<NetworKit::node> &getArticulationPointsOfComponent(NetworKit::count component) const {
        assert(has_run_);
        return articulation_points_of_component_[component];
    }

    bool isArticulationPoint(NetworKit::node u) const {
        assert(has_run_);
        return biconnected_components_.getComponentsOfNode(u).size() > 1;
    }

private:
    const NetworKit::Graph *const graph_;
    bool has_run_ = false;

    NetworKit::BiconnectedComponents biconnected_components_;

    std::vector<NetworKit::node> articulation_points_;
    std::vector<std::vector<NetworKit::node>> components_;
    std::vector<std::vector<NetworKit::node>> articulation_points_of_component_;
};

#endif // SMS_BLOCK_CUT_TREE_HPP
