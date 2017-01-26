
#ifndef TFTT_CELLREF_INT_H
#define TFTT_CELLREF_INT_H


#include "tftt.h"


namespace tftt {


//! Flags for a CellRef structure
enum CellRefFlags {
    TCR_COPYBOUNDARY = 1        //! If set, the cell is a copy boundary. Otherwise, the cell is a reflect boundary
};


//! A reference to a cell in the tree.
//! This is the public interface to most of the methods in the tree to be used by the program.
//!
//! A special case of the cell reference is if the group is null. In this case the reference is to a cell boundary,
//! the type of which is determined by the index member, using CellRefFlags (either copy or reflect).
struct CellRef {
public:
    //! Checks if the reference is to a boundary
    bool isBoundary() const;

    //! Checks if the reference is to a copy boundary
    bool isCopyBoundary() const;

    //! Checks if the reference is to a reflect boundary
    bool isReflectBoundary() const;

    // //! Reference to the parent cell
    // CellRef parent() const;

    // //! Returns the reference to the nth index child.
    // CellRef child(int n) const;
    // bool hasChildren() const;

    // //! The next cell in the space filling curve
    // CellRef next() const;

    // //! The previous cell in the space filling curve
    // CellRef prev() const;

    // //! The processor this cell will be evaluated on
    // int rank();

    // //! The identifier of the cell
    // ident_t id();

    // //! The data associated with the cell.
    // data_t& data();


    // void origin(double*& ret) const;
    // void centre(double*& ret) const;
    double width() const;
    double height() const;

};


typedef CellRef cell_t;


} // namespace tftt


#endif
