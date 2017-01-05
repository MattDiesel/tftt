

#ifndef TFTT_TREECELLREF_H
#define TFTT_TREECELLREF_H


#include "tftt/config.h"
#include "treegroup.h"


namespace tftt {


//! Flags for a TreeCellRef structure
enum TreeCellRefFlags {
    TCR_COPYBOUNDARY = 1        //! If set, the cell is a copy boundary. Otherwise, the cell is a reflect boundary
};


//! A reference to a cell in the tree.
//! This is the public interface to most of the methods in the tree to be used by the program.
//!
//! A special case of the cell reference is if the group is null. In this case the reference is to a cell boundary,
//! the type of which is determined by the index member, using TreeCellRefFlags (either copy or reflect).
struct TreeCellRef {
private:
    //! the group of siblings this cell belongs to
    TreeGroup* group;

    //! The index of the child within the sibling group.
    int index;

public:
    //! Ctor for a cell reference.
    //! Cell references should be generated from a call to Tree::find(), or by traversal of the tree.
    //!
    //! Note: With debug enabled, this constructor cannot be used to create boundary cells. Instead the
    //! TreeCellRef(bool) overload must be used.
    TreeCellRef(
            TreeGroup* gr,          //!< A pointer to the group
            int ind
        );

    //! Ctor for a reference to a boundary cell.
    TreeCellRef(
            bool copy               //!< If true, the cell is a copy boundary, otherwise the cell is a reflect boundary.
        );

    //! Checks if the reference is to a boundary
    bool isBoundary() const;

    //! Checks if the reference is to a copy boundary
    bool isCopyBoundary() const;

    //! Checks if the reference is to a reflect boundary
    bool isReflectBoundary() const;

    // //! Reference to the parent cell
    // TreeCellRef parent() const;

    // //! Returns the reference to the nth index child.
    // TreeCellRef child(int n) const;
    // bool hasChildren() const;

    // //! The next cell in the space filling curve
    // TreeCellRef next() const;

    // //! The previous cell in the space filling curve
    // TreeCellRef prev() const;

    // //! The processor this cell will be evaluated on
    // int rank();

    // //! The identifier of the cell
    // ident_t id();

    // //! The data associated with the cell.
    // data_t& data();


    // void origin(double*& ret) const;
    // void centre(double*& ret) const;
    // void size(double*& ret) const;

};


} // namespace tftt


#endif
