#include "mcp/auxiliary/statistics.hpp"

#include "networkit/clique/MaximalCliques.hpp"
#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/components/ConnectedComponents.hpp"

#include "mcp/auxiliary/small_ccs.hpp"

namespace mcp {

nlohmann::ordered_json graphMetricsJson(const NetworKit::Graph &g, bool fast) {
    nlohmann::ordered_json j;

    j["num nodes"] = g.numberOfNodes();
    j["num edges"] = g.numberOfEdges();

    // Degrees
    std::vector<unsigned int> degrees;
    for (auto v : g.nodeRange()) {
        degrees.push_back(g.degree(v));
    }

    auto degreeStats = IntStatistics(degrees);
    j.merge_patch(degreeStats.asJson());

    j["density"] =
        static_cast<double>(g.numberOfEdges()) * 2.0 / static_cast<double>(g.numberOfNodes() * (g.numberOfNodes() - 1));
    j["num degree 0"] = std::count(degrees.begin(), degrees.end(), 0);
    j["num degree 1"] = std::count(degrees.begin(), degrees.end(), 1);
    j["num degree 2"] = std::count(degrees.begin(), degrees.end(), 2);
    j["num degree 3"] = std::count(degrees.begin(), degrees.end(), 3);

    // Weights
    bool isInt = true;
    double intPart = 0;
    unsigned int num_zero = 0;
    std::vector<double> weights;
    double sumWeights = 0;
    double sumPositiveWeights = 0;
    for (auto e : g.edgeWeightRange()) {
        weights.push_back(e.weight);
        if (std::modf(e.weight, &intPart) != 0.0)
            isInt = false;
        if (e.weight == 0)
            num_zero++;
        if (e.weight > 0)
            sumPositiveWeights += e.weight;
        sumWeights += e.weight;
    }

    auto weightStats = FloatStatistics(weights);
    j.merge_patch(weightStats.asJson());

    j["weight type"] = isInt ? "integer" : "floating point";
    j["num zero weight edges"] = num_zero;
    j["sum of weights"] = sumWeights;
    j["sum of positive weights"] = sumPositiveWeights;

    // Components
    NetworKit::ConnectedComponents comps(g);
    comps.run();

    unsigned int num_one_comp = 0;
    std::vector<unsigned int> component_sizes;

    for (const auto &c : comps.getComponents()) {
        unsigned int size = c.size();

        if (size > 1)
            component_sizes.push_back(size);
        else
            num_one_comp++;
    }

    auto compStats = IntStatistics(component_sizes);

    j["number of cc"] = comps.getComponents().size();
    j["size of largest cc"] = compStats.max;
    j["size of smallest cc >1"] = compStats.min;
    j["mode of cc sizes"] = compStats.mode;
    j["number of size one cc"] = num_one_comp;

    NetworKit::BiconnectedComponents bi_cons(g);
    bi_cons.run();

    unsigned int num_two_bi_cons = 0;
    std::vector<unsigned int> bi_cons_sizes;

    for (const auto &c : bi_cons.getComponents()) {
        unsigned int size = c.size();

        if (size > 2)
            bi_cons_sizes.push_back(size);
        else if (size == 2)
            num_two_bi_cons++;
    }

    j["number of bcc"] = bi_cons.getComponents().size();

    if (not bi_cons_sizes.empty()) {
        auto BiConStats = IntStatistics(bi_cons_sizes);
        j["size of largest bcc"] = BiConStats.max;
        j["size of smallest bcc >2"] = BiConStats.min;
        j["mode of bcc sizes"] = BiConStats.mode;
        j["num size two bcc"] = num_two_bi_cons;
    }

    if (not fast) { // short chordless cycles
        SmallChordlessCycles ccs(g);
        ccs.run();
        j["num triangles"] = ccs.triangles.size();
        j["num 4-holes"] = ccs.fourHoles.size();

        auto mc = NetworKit::MaximalCliques(g);
        mc.run();
        auto cliques = mc.getCliques();

        std::vector<unsigned int> clique_sizes;

        for (const auto &i : cliques) {
            if (i.size() >= 5)
                clique_sizes.push_back(i.size());
        }

        if (!clique_sizes.empty()) {
            j["size biggest clique found of size >= 5"] = *std::max_element(clique_sizes.begin(), clique_sizes.end());
        }
    }

    return j;
}

void printGraphInformation(const NetworKit::Graph &g, std::ostream &out, bool fast) {
    out << "---------- Graph statistics -----------------" << std::endl;
    auto j = graphMetricsJson(g, fast);
    out << j.dump(4) << std::endl;
    out << "---------------------------------------------" << std::endl;
}

} // namespace mcp