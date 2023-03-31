#include <unordered_map>
#include <gtest/gtest.h>
#include <OpenVolumeMesh/Core/Handles.hh>

using namespace OpenVolumeMesh;

struct HandleTest : public ::testing::Test
{
};

TEST_F(HandleTest, unordered_map)
{
    auto m = std::unordered_map<VH, int>{};
    m[VH(0)] = 23;
    m[VH(1)] = 42;
    EXPECT_EQ(m[VH(0)], 23);
    EXPECT_EQ(m[VH(1)], 42);
}

