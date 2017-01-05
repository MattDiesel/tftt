
#include <gtest/gtest.h>

#include "tftt/config.h"
#include "tftt/treecellref.h"

using namespace tftt;


TEST(TreeCellRefTest, isBoundaryTest) {
    TreeCellRef ref(false); // reflect
    TreeCellRef cp(true); // copy
    TreeCellRef norm(reinterpret_cast<TreeGroup*>(0x1234), 3); // A "normal" reference

    // Check isBoundary
    ASSERT_EQ(ref.isBoundary(), true);
    ASSERT_EQ(cp.isBoundary(), true);
    ASSERT_EQ(norm.isBoundary(), false);

    // Check isCopyBoundary
    ASSERT_EQ(ref.isCopyBoundary(), false);
    ASSERT_EQ(cp.isCopyBoundary(), true);
    ASSERT_EQ(norm.isCopyBoundary(), false);

    // Check isReflectBoundary
    ASSERT_EQ(ref.isReflectBoundary(), true);
    ASSERT_EQ(cp.isReflectBoundary(), false);
    ASSERT_EQ(norm.isReflectBoundary(), false);
}

