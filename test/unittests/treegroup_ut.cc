
#include <gtest/gtest.h>

#include "tftt/tftt.h"
#include "tftt/tree.h"
#include "tftt/treegroup.h"



TEST(TreeGroupTest, CellIterator) {
    if (!tftt::gtree.root) {
        tftt::init(4.0, 2.0);
    }

    int ccells = 0;
    for (auto& c : *tftt::gtree.root) {
        ccells++;
    }
    ASSERT_EQ(ccells, 4);
}
