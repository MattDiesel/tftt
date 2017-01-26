

#ifndef TFTT_CELLREF_H
#define TFTT_CELLREF_H

#include "tftt.h"


namespace tftt {

struct TreeGroup;

//! Flags for a CellRef structure
enum CellRefFlags {
    TCR_INVALID = 0,
    TCR_COPYBOUNDARY = 1,        //! If set, the cell is a copy boundary. 
    TCR_REFLECTBOUNDARY = 2      //! If set, the cell is a reflect boundary
};


//! A reference to a cell in the tree.
//! This is the public interface to most of the methods in the tree to be used by the program.
//!
//! A special case of the cell reference is if the group is null. In this case the reference is to a cell boundary,
//! the type of which is determined by the index member, using CellRefFlags (either copy or reflect).
struct CellRef {
public:
    //! the group of siblings this cell belongs to
    TreeGroup* group;

    //! The index of the child within the sibling group.
    int index;

public:
    //! Ctor for a cell reference.
    //! Cell references should be generated from a call to Tree::find(), or by traversal of the tree.
    //!
    //! Note: With debug enabled, this constructor cannot be used to create boundary cells. Instead the
    //! CellRef(bool) overload must be used.
    CellRef(
            TreeGroup* gr,          //!< A pointer to the group
            int ind
        );

    //! Ctor for an invalid cell ref
    CellRef();

    //! Ctor for a reference to a boundary cell.
    CellRef(
            bool copy               //!< If true, the cell is a copy boundary, otherwise the cell is a reflect boundary.
        );

    //! Checks if the reference is valid.
    bool isValid() const;

    //! Checks if the reference is to a boundary
    bool isBoundary() const;

    //! Checks if the reference is to a copy boundary
    bool isCopyBoundary() const;

    //! Checks if the reference is to a reflect boundary
    bool isReflectBoundary() const;

    //! Reference to the parent cell
    CellRef parent() const;

    //! Returns the reference to the nth index child.
    CellRef child(int n) const;
    bool hasChildren() const;

    CellRef neighbour(int n) const;

    // //! The next cell in the space filling curve
    // CellRef next() const;

    // //! The previous cell in the space filling curve
    // CellRef prev() const;

    // //! The processor this cell will be evaluated on
    // int rank();

    // //! The identifier of the cell
    ident_t id() const;

    //! The data associated with the cell.
    data_t& data();


    double origin(int d) const;
    double centre(int d) const;
    double size(int d) const;

    bool operator==(const CellRef& rhs) const;
};


} // namespace tftt


#endif
