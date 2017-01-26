
#include <gtest/gtest.h>

#include "tftt/tftt.h"
#include "tftt/tree.h"


TEST(TfttTest, init) {
    tftt::init(4.0, 2.0);

    ASSERT_EQ(tftt::gtree.size[0], 4.0);
    ASSERT_EQ(tftt::gtree.size[1], 2.0);

    ASSERT_NE(tftt::gtree.root, nullptr);
    ASSERT_EQ(tftt::gtree.root->id, 0);
}


TEST(TfttTest, findTop) {
	if (!tftt::gtree.root) {
	    tftt::init(4.0, 2.0);
	}

	// Find cell {2}
    tftt::cell_t cl = tftt::find(2);

    ASSERT_EQ(cl.isValid(), true);
    ASSERT_EQ(cl.group, tftt::gtree.root);
    ASSERT_EQ(cl.index, 2);

	// Try find cell {5}
    cl = tftt::find(5);

    ASSERT_EQ(cl.isValid(), false);
}


TEST(TfttTest, cellRefBasics) {
    if (!tftt::gtree.root) {
        tftt::init(4.0, 2.0);
    }

    // Find cell {2}
    tftt::cell_t cl = tftt::find(2);

    ASSERT_EQ(cl.id().id, 2);
    ASSERT_EQ(cl.size(0), 2.0);
    ASSERT_EQ(cl.size(1), 1.0);
    ASSERT_EQ(cl.origin(0), 0.0);
    ASSERT_EQ(cl.origin(1), 1.0);
    ASSERT_EQ(cl.centre(0), 1.0);
    ASSERT_EQ(cl.centre(1), 1.5);
    ASSERT_EQ(cl.hasChildren(), false);

    tftt::cell_t parcl = cl.parent();

    ASSERT_EQ(parcl.isValid(), false);

    tftt::refine(cl);

    ASSERT_EQ(cl.hasChildren(), true);

    tftt::cell_t chc = tftt::find(cl.id().child(1));

    ASSERT_EQ(chc.isValid(), true);
    ASSERT_EQ(chc.size(0), 1.0);
    ASSERT_EQ(chc.size(1), 0.5);
    ASSERT_EQ(chc.origin(0), 1.0);
    ASSERT_EQ(chc.origin(1), 1.0);
    ASSERT_EQ(chc.hasChildren(), false);
}



TEST(TfttTest, cellRefFtt) {
    if (!tftt::gtree.root) {
        tftt::init(4.0, 2.0);
    }

    // Find cell {2}
    tftt::cell_t cl = tftt::find(2);

    ASSERT_EQ(cl.neighbour(0).isBoundary(), true);
    ASSERT_EQ(cl.neighbour(0).isCopyBoundary(), true);
    ASSERT_EQ(cl.neighbour(1).id(), 3);
    ASSERT_EQ(cl.neighbour(2).id(), 0);
    ASSERT_EQ(cl.neighbour(3).isBoundary(), true);
    ASSERT_EQ(cl.neighbour(3).isReflectBoundary(), true);
}
