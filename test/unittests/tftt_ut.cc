
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <gtest/gtest.h>

#include "tftt/tftt.h"
#include "tftt/tree.h"
#include "tftt/treegroup.h"

#include "util/formatstring.h"


using namespace util;


bool compareFiles(const std::string& p1, const std::string& p2) {
    std::ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
    std::ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }

    if (f1.tellg() != f2.tellg()) {
        return false; //size mismatch
    }

    //seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}



TEST(TfttTest, init) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    ASSERT_EQ(tftt::gtree.size[0], 4.0);
    ASSERT_EQ(tftt::gtree.size[1], 2.0);

    ASSERT_NE(tftt::gtree.root, nullptr);
    ASSERT_EQ(tftt::gtree.root->id, 0);

    ASSERT_EQ(tftt::gtree.root->neighbours[0].isBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[1].isBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[2].isBoundary(), true);
    ASSERT_EQ(tftt::gtree.root->neighbours[3].isBoundary(), true);

    ASSERT_EQ(tftt::gtree.root->neighbours[0].boundary(), 0);
    ASSERT_EQ(tftt::gtree.root->neighbours[1].boundary(), 1);
    ASSERT_EQ(tftt::gtree.root->neighbours[2].boundary(), 2);
    ASSERT_EQ(tftt::gtree.root->neighbours[3].boundary(), 3);
}


TEST(TfttTest, findTop) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    // Find cell {2}
    tftt::cell_t cl = tftt::find(2);

    ASSERT_EQ(cl.isValid(), true);
    ASSERT_EQ(cl.group, tftt::gtree.root);
    ASSERT_EQ(cl.index, 2);

    // Try find cell {5}
    cl = tftt::find(5);

    ASSERT_EQ(cl.isValid(), false);
}


TEST(TfttTest, atPosTop) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    double posValid[2] {1.0, 1.5};
    double posInValid[2] {5.0, 1.5};

    // Find cell {2}
    tftt::cell_t cl = tftt::atPos(posValid);

    ASSERT_EQ(cl.isValid(), true);
    ASSERT_EQ(cl.group, tftt::gtree.root);
    ASSERT_EQ(cl.index, 2);

    cl = tftt::atPos(posInValid);

    ASSERT_EQ(cl.isValid(), false);
}


TEST(TfttTest, atVertexTop) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    // Find cell {2}
    tftt::cell_t cl = tftt::atVertex(2);

    ASSERT_EQ(cl.isValid(), true);
    ASSERT_EQ(cl.group, tftt::gtree.root);
    ASSERT_EQ(cl.index, 2);
}


TEST(TfttTest, cellRefBasics) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    // Find cell {2}
    tftt::cell_t cl = tftt::find(2);

    ASSERT_EQ(cl.id().id, 2);
    ASSERT_EQ(cl.hasChildren(), false);
    ASSERT_EQ(cl.children(), nullptr);
    ASSERT_EQ(cl.level(), 0);

    tftt::cell_t parcl = cl.parent();

    ASSERT_EQ(parcl.isValid(), false);

    if (!cl.hasChildren())
        tftt::refine(cl);

    ASSERT_EQ(cl.children(), cl.child(0).group);
    ASSERT_EQ(cl.hasChildren(), true);

    tftt::cell_t chc = tftt::find(cl.id().child(1));

    ASSERT_EQ(chc.isValid(), true);
    ASSERT_EQ(chc.hasChildren(), false);
    ASSERT_EQ(chc.level(), 1);
}



TEST(TfttTest, cellRefFtt) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    // Find cell {2}
    tftt::cell_t cl = tftt::find(2);

    ASSERT_EQ(cl.neighbour(0).isBoundary(), true);
    ASSERT_EQ(cl.neighbour(1).id(), 3);
    ASSERT_EQ(cl.neighbour(2).id(), 0);
    ASSERT_EQ(cl.neighbour(3).isBoundary(), true);

    if (!cl.hasChildren())
        tftt::refine(cl);

    tftt::cell_t chc = cl.child(1);

    ASSERT_EQ(chc.neighbour(0), cl.child(0));
    ASSERT_EQ(chc.neighbour(3), cl.child(3));
    ASSERT_EQ(chc.neighbour(1), tftt::CellRef(tftt::gtree.root, 3));
    ASSERT_EQ(chc.neighbour(2), tftt::CellRef(tftt::gtree.root, 0));
}


TEST(TfttTest, drawMesh) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    tftt::cell_t cl = tftt::CellRef(tftt::gtree.root, 2);

    if (!cl.hasChildren())
        tftt::refine(cl);

    // Todo: Better test than just string comparison

    std::ostringstream oss;
    tftt::drawMesh(oss);

    std::string shouldBe = 
            "0 0\n2 0\n2 1\n0 1\n0 0\n\n"
            "2 0\n4 0\n4 1\n2 1\n2 0\n\n"
            "0 1\n1 1\n1 1.5\n0 1.5\n0 1\n\n"
            "1 1\n2 1\n2 1.5\n1 1.5\n1 1\n\n"
            "0 1.5\n1 1.5\n1 2\n0 2\n0 1.5\n\n"
            "1 1.5\n2 1.5\n2 2\n1 2\n1 1.5\n\n"
            "2 1\n4 1\n4 2\n2 2\n2 1\n\n";

    ASSERT_EQ(oss.str(), shouldBe);
}

TEST(TfttTest, TwoToOne) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    tftt::options.two2oneFlag = 1;

    tftt::cell_t cl = tftt::CellRef(tftt::gtree.root, 2);
    if (!cl.hasChildren()) {
        tftt::refine(cl);
    }

    tftt::cell_t chc = cl.child(1);
    if (!chc.hasChildren()) {
        tftt::refine(chc);
    }
    tftt::twoToOne(chc);

    // Should have refined 0 and 3 as well
    ASSERT_EQ(tftt::CellRef(tftt::gtree.root, 0).hasChildren(), true);
    ASSERT_EQ(tftt::CellRef(tftt::gtree.root, 1).hasChildren(), false);
    ASSERT_EQ(tftt::CellRef(tftt::gtree.root, 3).hasChildren(), true);

    ASSERT_EQ(tftt::CellRef(cl.children(), 0).hasChildren(), false);
    ASSERT_EQ(tftt::CellRef(cl.children(), 2).hasChildren(), false);
    ASSERT_EQ(tftt::CellRef(cl.children(), 3).hasChildren(), false);

}


TEST(TfttTest, insert) {
    tftt::reset();
    tftt::init(4.0, 2.0);

    tftt::ident_t idt = 0;
    for (int i = 0; i < 4; i++) {
        idt = idt.child(3);
    }

    // Find cell {2}
    tftt::options.two2oneFlag = 1;
    tftt::cell_t cl = tftt::insert(idt);

    ASSERT_EQ(cl.isValid(), true);

    tftt::drawMesh("testInsert.dat");
}



TEST(TfttTest, saveLoad) {
    std::string fileBefore = "testTree.before.dat";
    std::string fileAfter = "testTree.after.dat";

    tftt::reset();
    tftt::init(4.0, 2.0);

    tftt::ident_t idt = 0;
    for (int i = 0; i < 5; i++) {
        idt = idt.child(3);
    }

    tftt::options.two2oneFlag = 1;
    tftt::cell_t cl = tftt::insert(idt);
    ASSERT_EQ(cl.isValid(), true);

    tftt::saveTree("testTree.tr");
    tftt::drawMesh(fileBefore);
    tftt::reset();
    tftt::loadTree("testTree.tr");
    tftt::drawMesh(fileAfter);

    ASSERT_EQ(compareFiles(fileBefore, fileAfter), true);
}

TEST(TfttTest, splitToDisk) {
    std::string fileBefore = "testPartTree.before.dat";

    tftt::reset();
    tftt::init(4.0, 2.0);

    tftt::options.two2oneFlag = 2;

    tftt::ident_t idt = 0;
    for (int i = 0; i < 4; i++) {
        idt = idt.child(3);
    }

    tftt::cell_t cl = tftt::insert(idt);
    ASSERT_EQ(cl.isValid(), true);

    tftt::drawMesh(fileBefore);

    tftt::distribute(4);
    tftt::splitToDisk("testPartTree.r{0}.tr");

    for (int n = 0; n < 4; n++) {
        std::cout << "n = " << n << "\n";
        tftt::reset();
        tftt::loadTree(formatString("testPartTree.r{0}.tr", n), n);

        tftt::drawPartialMesh(formatString("testPartTree.r{0}.dat", n));
        tftt::drawGhosts(formatString("ghosts.r{0}.dat", n));

        std::cout << "\tGhosts: " << tftt::gtree.ghosts.size() << "\n";
    }

    // ASSERT_EQ(compareFiles(fileBefore, fileAfter), true);
}

TEST(TfttTest, interpFace) {
    tftt::fnData dtP = [](tftt::data_t& dt) {
        return dt.P;
    };
    double test;

    tftt::reset();
    tftt::init(1.0, 1.0);

    tftt::cell_t cl0 = tftt::find(0);
    tftt::cell_t cl1 = tftt::find(1);
    tftt::cell_t cl2 = tftt::find(2);
    tftt::cell_t cl3 = tftt::find(3);

    tftt::refine(cl0);
    tftt::refine(cl3);

    tftt::cell_t cl00 = cl0.child(0);
    tftt::refine(cl00);

    cl00.child(0)->P = 0;
    cl00.child(1)->P = 1;
    cl00.child(2)->P = 2;
    cl00.child(3)->P = 3;

    cl0.child(1)->P = 10;
    cl0.child(2)->P = 20;
    cl0.child(3)->P = 30;

    cl1->P = 100;
    cl2->P = 200;

    cl3.child(0)->P = 300;
    cl3.child(1)->P = 310;
    cl3.child(2)->P = 320;
    cl3.child(3)->P = 330;


    tftt::drawMesh("interpFace.mesh.dat");
    tftt::drawMatrix("interpFace.matrix.pgm", 400, 400, [](tftt::data_t& dt, int max) {
        return (dt.P / 400)*max;
    });

    // Case 1 - {3-2} to {3-3}
    test = tftt::interpFace(cl3.child(2), 1, dtP); 
    std::cout << "\t= " << test << "\n";

    // Case 2 - {2} to {3-*}
    test = tftt::interpFace(cl2, 1, dtP);
    std::cout << "\t= " << test << "\n";

    // Case 3 - {0-0-3} to {0-1}
    test = tftt::interpFace(cl00.child(3), 1, dtP);
    std::cout << "\t= " << test << "\n";

    // Case 4 - {0-3} to {1}
    test = tftt::interpFace(cl0.child(3), 1, dtP);
    std::cout << "\t= " << test << "\n";

    // Case 2 - {2} to {0-*}
    test = tftt::interpFace(cl2, 2, dtP);
    std::cout << "\t= " << test << "\n";

}
