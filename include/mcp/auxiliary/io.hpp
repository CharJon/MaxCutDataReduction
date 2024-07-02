#ifndef MCP_IO_HPP
#define MCP_IO_HPP

#include <filesystem>

#include "networkit/algebraic/MatrixTools.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/io/EdgeListReader.hpp"
#include "networkit/io/MatrixMarketReader.hpp"
#include "ogdf/fileformats/GraphIO.h"

namespace mcp {

/**
 * Interprets the @a matrix as adjacency matrix of a graph. If @a matrix is non-symmetric, the graph
 * will be directed.
 * @param matrix
 * @return The graph having an adjacency matrix equal to @a matrix.
 */
template <class Matrix>
NetworKit::Graph matrixToUndirectedGraphSafe(const Matrix &matrix) {
    NetworKit::Graph G(std::max(matrix.numberOfRows(), matrix.numberOfColumns()), true, false);
    matrix.forNonZeroElementsInRowOrder([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        if (weight != 1)
            throw std::runtime_error("Weighted graphs are not supported here.");
        G.addEdge(u, v, weight, true);
    });

    return G;
}

std::tuple<ogdf::Graph, std::unordered_map<NetworKit::node, ogdf::node>, std::unordered_map<int, NetworKit::node>>
nkToOgdf(const NetworKit::Graph &graph) {
    ogdf::Graph ogdfGraph;

    std::unordered_map<NetworKit::node, ogdf::node> nkToOgdfNodeMap;
    std::unordered_map<int, NetworKit::node> ogdfToNkNodeMap;

    for (auto u : graph.nodeRange()) {
        auto v = ogdfGraph.newNode();
        nkToOgdfNodeMap[u] = v;
        ogdfToNkNodeMap[v->index()] = u;
    }

    for (auto e : graph.edgeRange()) {
        auto u = e.u;
        auto v = e.v;
        ogdfGraph.newEdge(nkToOgdfNodeMap[u], nkToOgdfNodeMap[v]);
    }

    return {ogdfGraph, nkToOgdfNodeMap, ogdfToNkNodeMap};
}

NetworKit::Graph lazyGraphReader(const std::string &filePath, const std::string &format) {
    NetworKit::Graph g;
    if (format == "wel" || format == "wel0") {
        NetworKit::EdgeListReader r(' ', 0, "#", true, false);
        g = r.read(filePath);
    } else if (format == "wel1") {
        NetworKit::EdgeListReader r(' ', 1, "#", true, false);
        g = r.read(filePath);
    } else if (format == "mtx") {
        // read networkit graph from mtx file
        NetworKit::MatrixMarketReader mmr;
        auto res = mmr.read(filePath);
        g = matrixToUndirectedGraphSafe(res);
    } else if (format == "mc") {
        NetworKit::EdgeListReader r(' ', 1, "#", true, false);
        // copy content of mc file to a temporary .wel1 file
        auto tempDir = std::filesystem::temp_directory_path();
        auto tmpFile = tempDir / "tmp.wel1";
        std::ifstream src(filePath);
        std::ofstream dst(tmpFile);

        bool foundFirstNonCommentLine = false;
        std::string line;
        while (std::getline(src, line)) {
            if (foundFirstNonCommentLine) {
                dst << line << std::endl;
            } else {
                if (line[0] != '#') {
                    foundFirstNonCommentLine = true;
                } else {
                    dst << line << std::endl;
                }
            }
        }
        g = r.read(tmpFile.string());
        std::filesystem::remove(tmpFile);
    }

    return g;
}

void simpleMcIntWriter(const NetworKit::Graph &g, const std::string &path);

void simpleMcIntWriter(const NetworKit::Graph &g, const std::string &path) {
    std::ofstream file(path);

    if (!g.checkConsistency())
        throw std::runtime_error("Graph is inconsistent.");

    if (file.is_open()) {
        file << g.upperNodeIdBound() << " " << g.numberOfEdges() << std::endl;
        for (auto e : g.edgeWeightRange()) {
            assert(e.u <= e.v);

            if (e.weight == 0)
                throw std::runtime_error("Zero weight detected at edge: " + std::to_string(e.u) + " "
                                         + std::to_string(e.v));

            if (e.v == e.u)
                throw std::runtime_error("Selfloop detected at node: " + std::to_string(e.u));

            file << e.v + 1 << " " << e.u + 1 << " " << static_cast<int>(e.weight) << std::endl;
        }
    } else {
        throw std::runtime_error("File could not be opened: " + path + ".");
    }
    file.close();
}

} // namespace mcp

#endif // MCP_IO_HPP
