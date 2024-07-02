#include "tlx/container/d_ary_addressable_int_heap.hpp"
#include <gtest/gtest.h>

#include "mcp/auxiliary/active_elements.hpp"

template <class Type>
struct LessInVector {
    explicit LessInVector(const std::vector<Type> &vec) : vec_(&vec) {}
    bool operator()(uint64_t x, uint64_t y) const noexcept { return (*vec_)[x] < (*vec_)[y]; }

private:
    const std::vector<Type> *vec_;
};

TEST(activeElementsStack, empty) {
    auto stack = mcp::ActiveElementsStack(5u);
    EXPECT_TRUE(stack.empty());
}

TEST(activeElementsStack, fillAndEmpty) {
    auto stack = mcp::ActiveElementsStack(10u);
    stack.activate(1u);
    EXPECT_TRUE(stack.isActive(1u));
    stack.activate(3u);
    EXPECT_TRUE(stack.isActive(3u));
    ASSERT_FALSE(stack.empty());
    stack.deactivate(1u);
    ASSERT_FALSE(stack.empty());
    stack.deactivate(3u);
    EXPECT_TRUE(stack.empty());
}

TEST(activeElementsStack, tlx) {
    std::vector<uint8_t> vec = {1, 2, 3, 4, 5};
    auto s = tlx::DAryAddressableIntHeap<uint32_t, 2, LessInVector<uint8_t>>{LessInVector<uint8_t>{vec}};
    s.clear();
    s.push(1);
    s.push(2);
    s.push(3);
    ASSERT_TRUE(s.contains(1));
    ASSERT_TRUE(s.top() == 1);
    vec[1] = 10;
    s.update(1);
    ASSERT_TRUE(s.top() == 2);
}