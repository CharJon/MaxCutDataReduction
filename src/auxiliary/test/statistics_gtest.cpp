#include "networkit/graph/Graph.hpp"
#include <gtest/gtest.h>

#include "mcp/auxiliary/statistics.hpp"

TEST(statistics, simpleOdd1) {
    std::vector<int> values = {5, 2, 1};

    auto stats = mcp::IntStatistics(values);
    EXPECT_EQ(stats.median, 2);
}

TEST(statistics, simpleOdd2) {
    std::vector<int> values = {5, 2, 1, 3, 10};

    auto stats = mcp::IntStatistics(values);
    EXPECT_EQ(stats.median, 3);
}

TEST(statistics, simpleEven1) {
    std::vector<int> values = {5, 2};

    auto stats = mcp::IntStatistics(values);
    EXPECT_EQ(stats.median, 2);
}

TEST(statistics, simpleEven2) {
    std::vector<int> values = {5, 10, 6, 4, 3, 2, 6, 7, 9, 3};

    auto stats = mcp::IntStatistics(values);
    EXPECT_EQ(stats.median, 5);
}
