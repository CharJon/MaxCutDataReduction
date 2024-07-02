#ifndef SMS_CL_PARSER_HPP
#define SMS_CL_PARSER_HPP

#include "cxxopts.hpp"
#include <optional>

namespace mcp {

class CLParser {
public:
    CLParser(int argc, char **argv);

    std::string getFileName() const;

    std::string getType() const;

    int getSeed() const;

    double getTimeLimit() const;

    std::optional<std::string> getStatsFile() const;

    int getThreads() const;

    int getConnectivity() const;

    std::optional<std::string> getOutputFile() const;

    unsigned int getTriangleFlag() const;

    int getCliqueLevel() const;

    bool getStableReduceFlag() const;

    bool getSolve() const;

    bool getHeuristic() const;

    std::optional<int> getScalingFactor() const;

private:
    cxxopts::ParseResult parseResult_;

    static cxxopts::Options createParser();
};

} // namespace mcp

#endif // SMS_CL_PARSER_HPP
