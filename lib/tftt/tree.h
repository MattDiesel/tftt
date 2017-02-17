
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

    bool destroying;
    ~Tree();

	// For thread
    cell_t first;
    cell_t last;
    cell_t firstActive;
    cell_t lastActive;

	//! The set of ghost groups on this processor
	// std::set<TreeGroup*> ghosts;

    node_t rank;

};


//! The tree instance
extern Tree gtree;


} // namespace tftt


#endif
