
#include <gtest/gtest.h>

#include "tftt/tftt.h"
#include "tftt/tree.h"
#include "tftt/treegroup.h"


TEST(TfttTest, init) {
    if (!tftt::gtree.root) {
        tftt::init(4.0, 2.0);
    }

    ASSERT_EQ(tftt::gtree.size[0], 4.0);
    ASSERT_EQ(tftt::gtree.size[1], 2.0);

    ASSERT_NE(tftt::gtree.root, nullptr);
    ASSERT_EQ(tftt::gtree.root->id, 0);

    ASSERT_EQ(tftt::gtree.root->neighbours[0].isBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[1].isBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[2].isBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[3].isBoundary(), true);

    ASSERT_EQ(tftt::gtree.root->neighbours[0].isCopyBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[1].isCopyBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[0].isReflectBoundary(), false);
    ASSERT_EQ(tftt::gtree.root->neighbours[1].isReflectBoundary(), false);

    ASSERT_EQ(tftt::gtree.root->neighbours[2].isReflectBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[3].isReflectBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[2].isCopyBoundary(), false);
    ASSERT_EQ(tftt::gtree.root->neighbours[3].isCopyBoundary(), false);
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
    ASSERT_EQ(cl.children(), nullptr);

    tftt::cell_t parcl = cl.parent();

    ASSERT_EQ(parcl.isValid(), false);

    tftt::refine(cl);

    ASSERT_EQ(cl.children(), cl.child(0).group);
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

    if (!cl.hasChildren())
        tftt::refine(cl);

    tftt::cell_t chc = cl.child(1);

    ASSERT_EQ(chc.neighbour(0), cl.child(0));
    ASSERT_EQ(chc.neighbour(3), cl.child(3));
    ASSERT_EQ(chc.neighbour(1), tftt::CellRef(tftt::gtree.root, 3));
    ASSERT_EQ(chc.neighbour(2), tftt::CellRef(tftt::gtree.root, 0));
}
