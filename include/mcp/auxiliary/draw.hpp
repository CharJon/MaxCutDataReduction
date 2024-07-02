#ifndef MCP_DRAW_HPP
#define MCP_DRAW_HPP

#include "ogdf/energybased/FMMMLayout.h"
#include "ogdf/basic/GraphAttributes.h"
#include "ogdf/fileformats/GraphIO.h"

void drawSvg(const ogdf::Graph &g, const std::string &path) {
    ogdf::GraphAttributes GA2(g);
    for (ogdf::node v: g.nodes)
        GA2.width(v) = GA2.height(v) = 5.0;

    ogdf::FMMMLayout fmmm;

    fmmm.useHighLevelOptions(true);
    fmmm.unitEdgeLength(15.0);
    fmmm.newInitialPlacement(true);
    fmmm.qualityVersusSpeed(ogdf::FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);

    fmmm.call(GA2);
    ogdf::GraphIO::write(GA2, path, ogdf::GraphIO::drawSVG);
}

#endif //MCP_DRAW_HPP
