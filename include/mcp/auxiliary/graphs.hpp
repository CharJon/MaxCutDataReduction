#ifndef SMS_GRAPHS_HPP
#define SMS_GRAPHS_HPP

#include <concepts>
#include <vector>

#include "networkit/components/BiconnectedComponents.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"

namespace mcp {

using NetworKit::count;
using NetworKit::Edge;
using NetworKit::edgeweight;
using NetworKit::Graph;
using NetworKit::node;

std::tuple<node, count> maxDegreeNode(const NetworKit::Graph &g);

unsigned int numDegreeNodes(const NetworKit::Graph &graph, unsigned int degree);

unsigned int numWeightEdges(const Graph &graph, edgeweight weight);

Graph unweightedComplementGraph(const Graph &graph);

void scaleEdgeWeights(Graph &graph, double factor);

} // namespace mcp

#endif // SMS_GRAPHS_HPP
