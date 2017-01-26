
#include <sstream>

#include <gtest/gtest.h>

#define DIM 2

#include "tftt/tftt.h"
#include "tftt/treeid.h"


using namespace tftt;

typedef TreeId<uint64_t> qtidl_t;


TEST(TreeIdTest, BasicTests) {
	qtidl_t myId;
	ASSERT_EQ(myId.level(), 0);

	qtidl_t myChild = myId.child(1);
	ASSERT_EQ(myChild.level(), 1);
	ASSERT_EQ(myChild.orthant(), 1);
	ASSERT_EQ(myChild.parent(), myId);

	qtidl_t myGrandChild = myChild.child(2);
	ASSERT_EQ(myGrandChild.level(), 2);
	ASSERT_EQ(myGrandChild.orthant(), 2);
	ASSERT_EQ(myGrandChild.parent(), myChild);
}


TEST(TreeIdTest, PathPrinting) {
	std::ostringstream oss;

	qtidl_t myId;
	myId = myId.child(1);
	myId = myId.child(2);
	myId = myId.child(3);
	myId = myId.child(0);
	myId = myId.child(1);
	myId = myId.child(2);
	myId = myId.child(3);

	oss << myId;

	ASSERT_EQ(oss.str(), "0-1-2-3-0-1-2-3");
}


// TEST(TreeIdTest, BuildConfig) {
// 	tftt::DIM = 2;
// 	tftt::ident_t myId;

// 	int cbytes = sizeof(myId.id);
// 	ASSERT_GE((cbytes-1)*8 / tftt::DIM, tftt::MaxDepth) << "ID Datatype insufficient to store maximum depth of the tree.";
// }


// TEST(TreeIdTest, BuildSpecific) {
// 	tftt::ident_t myId;
// 	ASSERT_EQ(myId.level(), 0);

// 	tftt::ident_t myChild = myId.child(1);
// 	ASSERT_EQ(myChild.level(), 1);
// 	ASSERT_EQ(myChild.orthant(), 1);
// 	ASSERT_EQ(myChild.parent(), myId);

// 	tftt::ident_t myGrandChild = myChild.child(2);
// 	ASSERT_EQ(myGrandChild.level(), 2);
// 	ASSERT_EQ(myGrandChild.orthant(), 2);
// 	ASSERT_EQ(myGrandChild.parent(), myChild);

// 	std::ostringstream oss;

// 	myId = myId.child(1);
// 	myId = myId.child(2);
// 	myId = myId.child(3);

// 	oss << myId;

// 	ASSERT_EQ(oss.str(), "0-1-2-3");
// }
