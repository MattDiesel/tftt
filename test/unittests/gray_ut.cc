
#include <gtest/gtest.h>

#include "tftt/gray.h"

using namespace tftt::utils;

TEST(GrayTest, basic) {
    ASSERT_EQ(toGray(0), 0);
    ASSERT_EQ(toGray(1), 1);
    ASSERT_EQ(toGray(2), 3);
    ASSERT_EQ(toGray(3), 2);
    ASSERT_EQ(toGray(7), 4);
}
