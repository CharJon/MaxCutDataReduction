#ifndef MCP_ACTIVE_ELEMENTS_HPP
#define MCP_ACTIVE_ELEMENTS_HPP

#include <cassert>
#include <concepts>
#include <limits>
#include <vector>

namespace mcp {

template <std::unsigned_integral T>
class ActiveElementsStack {
public:
    explicit ActiveElementsStack(T n) : n_(n), stack_(), posInStack_(n_, kINACTIVE) { stack_.reserve(n_); }

    bool isActive(T x) {
        assert(x < n_);
        return posInStack_[x] != kINACTIVE;
    }

    void activate(T x) {
        assert(x < n_);
        if (!isActive(x)) {
            posInStack_[x] = stack_.size();
            stack_.push_back(x);
        }
    }

    void deactivate(T x) {
        assert(x < n_);
        if (isActive(x)) {
            auto back = stack_.back();
            posInStack_[back] = posInStack_[x];
            stack_[posInStack_[back]] = back;
            stack_.pop_back();
            posInStack_[x] = kINACTIVE;
        }
    }

    T popBack() {
        assert(!empty());
        auto last = stack_.back();
        posInStack_[last] = kINACTIVE;
        stack_.pop_back();
        return last;
    }

    bool empty() { return stack_.empty(); }

private:
    constexpr static T kINACTIVE = std::numeric_limits<T>::max();

    T n_;
    std::vector<T> stack_;
    std::vector<T> posInStack_;
};

} // namespace mcp

#endif // MCP_ACTIVE_ELEMENTS_HPP
