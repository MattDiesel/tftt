

#ifndef TFTT_CELLREF_H
#define TFTT_CELLREF_H


#include <ostream>

#include "tftt.h"


namespace tftt {

struct TreeGroup;

//! Flags for a CellRef structure
enum CellRefFlags {
    TCR_INVALID = 0,
    TCR_ROOT = 1
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

    //! Ctor for a cell ref with flags
    //! e.g. CellRef(1) gives the root cell.
    CellRef(int flag);

    //! Checks if the reference is valid.
    bool isValid() const;

    //! Checks if the reference is valid.
    bool isRoot() const;

    //! Checks if the reference is to a boundary
    bool isBoundary() const;
    int boundary() const;

    //! Reference to the parent cell
    CellRef parent() const;

    //! Returns the reference to the nth index child.
    CellRef child(int n) const;
    CellRef childOnFace(int fc, int n) const;
    bool hasChildren() const;
    bool hasGrandChildren() const;
    TreeGroup* children() const;
    // TreeGroup* group() const;
    // TreeBoundaryGroup* bgroup() const;

    CellRef neighbour(int n) const;

    //! The orientation of the space filling curve at this cell. 
    int orientation() const;

    //! The next cell in the space filling curve
    CellRef next() const;

    //! The previous cell in the space filling curve
    CellRef prev() const;

    bool isLastInGroup() const;
    bool isFirstInGroup() const;

    //! The processor this cell will be evaluated on
    node_t& rank();
    node_t const& rank() const;

    // //! The identifier of the cell
    ident_t id() const;

    int level() const;

    //! The data associated with the cell.
    data_t& data();
    data_t const& data() const;

    //! Average values of all children (recursive)
    double avrChildren(fnData dt) const;

    //! Interpolate a "neighbours" data
    double ngbVal(int nb, fnData dt, double* ifBoundary = nullptr) const;

    //! The data associated with a cell face. 
    facedata_t& facedata(int dir);
    facedata_t const& facedata(int dir) const;


    double origin(int d) const;
    double centre(int d) const;
    double size(int d) const;
    double vertex(int v, int d) const;
    bool containsPoint(double pt[DIM]) const;

    void sizes(double ret[DIM]) const;
    void origins(double ret[DIM]) const;
    void vertices(double ret[1<<DIM][DIM]) const;

    bool operator==(const CellRef& rhs) const;
    bool operator!=(const CellRef& rhs) const;
    bool operator<(const CellRef& rhs) const;

    data_t* operator->();
    data_t const* operator->() const;

    std::string path() const;
};


} // namespace tftt


std::ostream& operator<<(std::ostream& os, const tftt::CellRef& cr);


#endif
