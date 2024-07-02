#include "mcp/auxiliary/graphs.hpp"

#include <map>
#include <vector>

#include "networkit/graph/Graph.hpp"

namespace mcp {

using NetworKit::Graph;
using NetworKit::node;

std::tuple<node, count> maxDegreeNode(const NetworKit::Graph &g) {
    NetworKit::node maxDegreeNode = 0;
    NetworKit::count maxDegree = g.degree(maxDegreeNode);
    for (auto u : g.nodeRange()) {
        NetworKit::count currentDegree = g.degree(u);
        if (currentDegree > maxDegree) {
            maxDegree = currentDegree;
            maxDegreeNode = u;
        }
    }
    return {maxDegreeNode, maxDegree};
}

unsigned int numDegreeNodes(const Graph &graph, unsigned int degree) {
    auto nr = graph.nodeRange();
    return std::count_if(nr.begin(), nr.end(),
                         [&graph, degree](NetworKit::node v) { return graph.degree(v) == degree; });
}

unsigned int numWeightEdges(const Graph &graph, edgeweight weight) {
    auto er = graph.edgeWeightRange();
    return std::count_if(er.begin(), er.end(), [weight](NetworKit::WeightedEdge e) { return e.weight == weight; });
}

Graph unweightedComplementGraph(const Graph &graph) {
    assert(!graph.isDirected());
    assert(graph.numberOfNodes() == graph.upperNodeIdBound());

    auto complementGraph = Graph(graph.upperNodeIdBound());

    auto lastNeighbor = std::vector<NetworKit::node>(
        graph.upperNodeIdBound(),
        std::numeric_limits<NetworKit::node>::max()); // max of unsigned type here allows us to get 0 if we add 1
    for (NetworKit::node u = 0; u < graph.upperNodeIdBound(); u++) {
        for (auto neighbor : graph.neighborRange(u)) {
            for (auto v = lastNeighbor[neighbor] + 1; v < std::min(u, neighbor); v++) {
                assert(!complementGraph.hasEdge(neighbor, v));
                complementGraph.addEdge(neighbor, v);
            }
            lastNeighbor[neighbor] = u;
        }
    }

    for (NetworKit::node u = 1; u < graph.upperNodeIdBound(); u++) {
        for (NetworKit::node v = lastNeighbor[u] + 1; v < u; v++) {
            complementGraph.addEdge(u, v);
        }
    }

    return complementGraph;
}

void scaleEdgeWeights(Graph &graph, double factor) {
    for (auto e : graph.edgeWeightRange()) {
        auto newWeight = round(e.weight * factor);
        graph.setWeight(e.u, e.v, newWeight);
    }
}

} // namespace mcp
