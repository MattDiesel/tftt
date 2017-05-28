
#ifndef TFTT_TREEGROUP_H
#define TFTT_TREEGROUP_H

#include <iterator>
#include <array>

#include "../config.h"
#include "../treeid.h"
#include "../cellref.h"
#include "treecell.h"


namespace tftt {


//! A group of adjacent sibling cells. (A quad in 2d, oct in 3d)
//! This reduces the memory requirement by not having to link all cells to their neighbours.
//!
//! TreeGroups should never be directly edited by the user.
struct TreeGroup {

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
    TreeGroup(CellRef parent);
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
    std::array<CellRef,DIM*2> neighbours;

    //! Pointer to the parent group.
    CellRef parent;

    // For Thread

    //! The Hilbert orientation of this group
    int orientation;

    //! The next group in the space filling curve.
    CellRef next;

    //! The previous group in the space filling curve
    CellRef prev;


    class cell_iterator {
    public:
        CellRef cr;

        cell_iterator(TreeGroup* gr, int ind)
            : cr(gr, ind) {
        }

        cell_iterator operator++() {
            cell_iterator i = *this;
            cr.index++;
            return i;
        }
        cell_iterator operator++(int junk) {
            cr.index++;
            return *this;
        }
        CellRef& operator*() {
            return cr;
        }
        CellRef* operator->() {
            return &cr;
        }
        bool operator==(const cell_iterator& rhs) {
            return cr == rhs.cr;
        }
        bool operator!=(const cell_iterator& rhs) {
            return !(cr == rhs.cr);
        }
    };

    cell_iterator begin() {
        return cell_iterator(this, 0);
    }
    cell_iterator end() {
        return cell_iterator(this, 1<<DIM);
    }
};


} // namespace tftt


#endif
