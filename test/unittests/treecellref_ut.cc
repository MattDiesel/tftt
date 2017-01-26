
#include <gtest/gtest.h>

#include "tftt/tftt.h"
#include "tftt/cellref.h"
#include "tftt/tree.h"


using namespace tftt;


TEST(TreeCellRefTest, isBoundaryTest) {
    CellRef ref(false); // reflect
    CellRef cp(true); // copy
    CellRef norm(reinterpret_cast<TreeGroup*>(0x1234), 3); // A "normal" reference

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


TEST(TreeCellRefTest, locations) {
    if (!tftt::gtree.root) {
        tftt::init(4.0, 2.0);
    }

    // Find cell {2}
    tftt::cell_t cl = tftt::find(2);

    ASSERT_EQ(cl.origin(0), 0.0);
    ASSERT_EQ(cl.origin(1), 1.0);
    ASSERT_EQ(cl.centre(0), 1.0);
    ASSERT_EQ(cl.centre(1), 1.5);
    ASSERT_EQ(cl.size(0), 2.0);
    ASSERT_EQ(cl.size(1), 1.0);
    ASSERT_EQ(cl.vertex(0, 0), 0.0);
    ASSERT_EQ(cl.vertex(0, 1), 1.0);
    ASSERT_EQ(cl.vertex(1, 0), 2.0);
    ASSERT_EQ(cl.vertex(1, 1), 1.0);
    ASSERT_EQ(cl.vertex(2, 0), 0.0);
    ASSERT_EQ(cl.vertex(2, 1), 2.0);
    ASSERT_EQ(cl.vertex(3, 0), 2.0);
    ASSERT_EQ(cl.vertex(3, 1), 2.0);
}
