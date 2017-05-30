
#ifndef TFTT_TREEGROUP_H
#define TFTT_TREEGROUP_H

#include <iterator>
#include <array>

#include "../config.h"
#include "../treeid.h"
#include "treecell.h"


namespace tftt {


//! A group of adjacent sibling cells. (A quad in 2d, oct in 3d)
//! This reduces the memory requirement by not having to link all cells to their neighbours.
//!
//! TreeGroups should never be directly edited by the user.
struct TreeGroup {
    class cell_iterator;

    bool isBoundary() const {
        return id.isBoundary();
    };
    int boundary() const {
        return id.boundary();
    }
    void setBoundary(int b) {
        id = ident_t::boundary(b);
    }

    // Base Tree
    TreeGroup();
    TreeGroup(int b);
    TreeGroup(cell_t parent);
    ~TreeGroup();

    //! The identifier of the group, equal to the identifier of the first cell.
    //! The id of other cells in the group can be found by adding their index to this id.
    ident_t id;

    //! The origin (top left corner in 2d) of the group.
    double origin[DIM];

    //! The cells in this group
    std::array<TreeCell, 1<<DIM> cells;

    // FTT

    //! Pointers to the neighbouring groups
    std::array<cell_t,DIM*2> neighbours;

    //! Pointer to the parent group.
    cell_t parent;

    // For Thread

    //! The Hilbert orientation of this group
    int orientation;

    //! The next group in the space filling curve.
    cell_t next;

    //! The previous group in the space filling curve
    cell_t prev;

    cell_iterator begin();
    cell_iterator end();


    class cell_iterator {
    public:
        cell_t cr;

        cell_iterator(TreeGroup* gr, int ind);
        cell_iterator operator++();
        cell_iterator operator++(int junk);
        cell_t& operator*();
        cell_t* operator->();
        bool operator==(const cell_iterator& rhs);
        bool operator!=(const cell_iterator& rhs);
    };
};


} // namespace tftt


#endif
