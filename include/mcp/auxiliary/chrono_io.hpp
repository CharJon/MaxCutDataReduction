#ifndef MCP_CHRONO_IO_HPP
#define MCP_CHRONO_IO_HPP

#include <chrono>
#include <concepts>
#include <iostream>
#include <type_traits>

template <typename T>
inline std::string suffix()
    requires std::same_as<T, std::milli>
{
    return "ms";
}

template <typename T>
inline std::string suffix()
    requires std::same_as<T, std::micro>
{
    return "µs";
}

template <typename T>
inline std::string suffix()
    requires std::same_as<T, std::nano>
{
    return "ns";
}

template <typename Rep>
std::ostream &operator<<(std::ostream &os, const std::chrono::duration<Rep> &duration) {
    os << duration.count() << "s";
    return os;
}

template <typename Rep, typename Ratio>
std::ostream &operator<<(std::ostream &os, const std::chrono::duration<Rep, Ratio> &duration) {
    os << duration.count() << suffix<Ratio>();
    return os;
}

#endif // MCP_CHRONO_IO_HPP
