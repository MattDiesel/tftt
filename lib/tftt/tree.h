
#ifndef TFTT_TREE_H
#define TFTT_TREE_H

#include <set>
#include <vector>

#include "tftt.h"
#include "treegroup.h"


namespace tftt {


struct Tree {
    TreeGroup* root;
    size_t ccells;
    size_t cactive;

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
    std::vector<std::set<cell_t, crparless>> ghosts;
    std::vector<std::set<cell_t, crparless>> borders;

    std::vector<std::vector<cell_t>> rawGhosts;
    std::vector<std::vector<cell_t>> rawBorders;

    std::vector<data_t*> ghostData;
    std::vector<data_t*> borderData;

    node_t rank;

};


//! The tree instance
extern Tree gtree;


} // namespace tftt


#endif
