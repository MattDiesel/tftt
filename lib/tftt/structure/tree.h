
#ifndef TFTT_TREE_H
#define TFTT_TREE_H


#include <set>
#include <vector>

#include "../config.h"
#include "../cellref.h"
#include "treegroup.h"


namespace tftt {


struct Tree {
    TreeGroup* root;
    size_t ccells;
    size_t cactive;

    double size[DIM];

    TreeGroup* boundGroups;

    bool isNeuman[2*DIM];
    double dirichletValue[2*DIM];

    bool destroying;
    ~Tree();

    // For thread
    cell_t first;
    cell_t last;
    cell_t firstActive;
    cell_t lastActive;

    //! The set of ghost groups on this processor
    std::vector<std::set<cell_t, cell_t::parless>> ghosts;
    std::vector<std::set<cell_t, cell_t::parless>> borders;

    std::vector<std::vector<cell_t>> rawGhosts;
    std::vector<std::vector<cell_t>> rawBorders;

    std::vector<data_t*> ghostData;
    std::vector<data_t*> borderData;

    std::vector<uint32_t*> ghostAdaptVectors;
    std::vector<uint32_t*> borderAdaptVectors;

    node_t rank;
};


//! The tree instance
extern Tree gtree;


} // namespace tftt


#endif
