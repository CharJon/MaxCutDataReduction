#include "mcp/auxiliary/cl_parser.hpp"

#include "cxxopts.hpp"
#include <optional>

namespace mcp {
#define INSTANCE "instance"
#define TYPE "type"
#define SEED "seed"
#define TIMELIMIT "timelimit"
#define STATSFILE "statsfile"
#define OUTPUTFILE "outputfile"
#define THREADS "threads"
#define CONNECTIVITY "connectivity"
#define TRIANGLES "triangles"
#define CLIQUELEVEL "clique-level"
#define STABLER "stable-reduce"
#define SCALE "scale"
#define SOLVE "solve"

cxxopts::Options CLParser::createParser() {
    cxxopts::Options options("SMS", "Max cut solver based on SCIP  ");

    options.add_options()
        // positional
        (INSTANCE, "Instance file name", cxxopts::value<std::string>())(TYPE, "Type (optional)",
                                                                        cxxopts::value<std::string>())
        // pseudo-optional (have default values)
        (SEED, "Random seed", cxxopts::value<int>()->default_value("-1"))                       //
        (TIMELIMIT, "Maximum time in seconds", cxxopts::value<double>()->default_value("1e20")) //
        (THREADS, "Number of threads", cxxopts::value<int>()->default_value("1"))               //
        (CONNECTIVITY, "Up to which connectivity shall the preprocessing try to contract subgraphs",
         cxxopts::value<int>()->default_value("3"))                                                           //
        (CLIQUELEVEL, "Clique level to use", cxxopts::value<int>()->default_value("2"))                       //
        (TRIANGLES, "Which triangle rules shall be used", cxxopts::value<unsigned int>()->default_value("7")) //
        (STABLER, "Shall stable reduce be used", cxxopts::value<bool>()->default_value("true"))               //
        (SOLVE, "Which solver to use for remaining graphs (exact, heuristic, none)",
         cxxopts::value<std::string>()->default_value("false")) //
        // real optional
        (STATSFILE, "File to store compact stats in", cxxopts::value<std::string>()) //
        (OUTPUTFILE, "File to store output in", cxxopts::value<std::string>())       //
        (SCALE, "Multiply and round", cxxopts::value<int>())                         //
        ;
    options.parse_positional({INSTANCE, TYPE});

    return options;
}

CLParser::CLParser(int argc, char **argv) {
    parseResult_ = createParser().parse(argc, argv);
    if (!(parseResult_[SOLVE].as<std::string>() == "exact" || parseResult_[SOLVE].as<std::string>() == "heuristic"
          || parseResult_[SOLVE].as<std::string>() == "none")) {
        throw std::runtime_error("Solve flag unknown, was " + parseResult_[SOLVE].as<std::string>() + ".");
    }
}

std::string CLParser::getFileName() const {
    if (parseResult_.count(INSTANCE) == 0) {
        return {};
    } else {
        return parseResult_[INSTANCE].as<std::string>();
    }
}

std::string CLParser::getType() const {
    if (parseResult_.count(TYPE) == 0) {
        uint64_t positionLastDot = getFileName().find_last_of('.');
        return getFileName().substr(positionLastDot + 1);
    } else {
        return parseResult_[TYPE].as<std::string>();
    }
}

int CLParser::getSeed() const {
    return parseResult_[SEED].as<int>();
}

double CLParser::getTimeLimit() const {
    return parseResult_[TIMELIMIT].as<double>();
}

std::optional<std::string> CLParser::getStatsFile() const {
    if (parseResult_.count(STATSFILE) == 0)
        return std::nullopt;
    else
        return parseResult_[STATSFILE].as<std::string>();
}

std::optional<std::string> CLParser::getOutputFile() const {
    if (parseResult_.count(OUTPUTFILE) == 0)
        return std::nullopt;
    else
        return parseResult_[OUTPUTFILE].as<std::string>();
}

int CLParser::getThreads() const {
    return parseResult_[THREADS].as<int>();
}

int CLParser::getConnectivity() const {
    return parseResult_[CONNECTIVITY].as<int>();
}

unsigned int CLParser::getTriangleFlag() const {
    return parseResult_[TRIANGLES].as<unsigned int>();
}

int CLParser::getCliqueLevel() const {
    return parseResult_[CLIQUELEVEL].as<int>();
}

bool CLParser::getStableReduceFlag() const {
    return parseResult_[STABLER].as<bool>();
}

bool CLParser::getSolve() const {
    return parseResult_[SOLVE].as<std::string>() == "exact";
}

bool CLParser::getHeuristic() const {
    return parseResult_[SOLVE].as<std::string>() == "heuristic";
}

std::optional<int> CLParser::getScalingFactor() const {
    if (parseResult_.count(SCALE) == 0)
        return {};
    else
        return parseResult_[SCALE].as<int>();
}

} // namespace mcp
