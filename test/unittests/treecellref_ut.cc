
#include <gtest/gtest.h>
#include <cmath>

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

TEST(TreeCellRefTest, dataAccess) {
    if (!tftt::gtree.root) {
        tftt::init(4.0, 2.0);
    }

    tftt::cell_t cl = tftt::find(2);

    cl.data().P = 1337;

    ASSERT_EQ(cl.data().P, 1337);
    ASSERT_EQ(cl->P, 1337);

    tftt::cell_t const& cl_const = cl;

    ASSERT_EQ(cl_const.data().P, 1337);
    ASSERT_EQ(cl_const->P, 1337);

    // Face data

    cl.facedata(0).poisCoef = 1337;
    ASSERT_EQ(cl.facedata(0).poisCoef, 1337);
}

TEST(TreeCellRefTest, interpChild) {

    tftt::reset();
    tftt::init(4.0, 2.0);

    tftt::cell_t cl0 = tftt::find(0);
    tftt::cell_t cl1 = tftt::find(1);
    tftt::cell_t cl2 = tftt::find(2);
    tftt::cell_t cl3 = tftt::find(3);

    tftt::refine(cl0);

    cl0.child(0)->P = 0.0;
    cl0.child(1)->P = 0.1;
    cl0.child(2)->P = 0.2;
    cl0.child(3)->P = 0.3;

    cl1->P = 1.0;
    cl2->P = 2.0;
    cl3->P = 3.0;

    double c12 = tftt::interpChild(cl1, 2, 0, [](tftt::data_t& dt) {
        return dt.P;
    });

    ASSERT_LE(std::abs(c12 - 1.1), 1e-4);


    tftt::refine(cl3);

    cl3.child(0)->P = 3.0;
    cl3.child(1)->P = 3.1;
    cl3.child(2)->P = 3.2;
    cl3.child(3)->P = 3.3;

    tftt::drawMatrix("interpChild.matrix.pgm", 800, 400, [](tftt::data_t& dt, int max) {
        return (dt.P / 4.0)*max;
    });

    c12 = tftt::interpChild(cl1, 2, 0, [](tftt::data_t& dt) {
        return dt.P;
    });

    ASSERT_LE(std::abs(c12 - 1.22222), 1e-4);
}
