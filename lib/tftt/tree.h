
#ifndef TFTT_TREE_H
#define TFTT_TREE_H

#include <set>

#include "tftt.h"
#include "treegroup.h"


namespace tftt {


struct Tree {
	TreeGroup* root;
	size_t ccells;

	double size[DIM];

    TreeGroup* boundGroups;

    ~Tree();

	// For thread
	// TreeGroup* first;
	// TreeGroup* last;

	//! The set of ghost groups on this processor
	// std::set<TreeGroup*> ghosts;

};


//! The tree instance
extern Tree gtree;


} // namespace tftt


#endif
