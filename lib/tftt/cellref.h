

#ifndef TFTT_CELLREF_H
#define TFTT_CELLREF_H


#include <string>
#include <ostream>

#include "config.h"
#include "treeid.h"


namespace tftt {


struct TreeCell;
struct TreeGroup;
struct tagNeighbours;
struct tagPoissonNeighbours;


//! Flags for a CellRef structure
enum CellRefFlags {
    TCR_INVALID = 0,
    TCR_ROOT = 1
};


//! A reference to a cell in the tree.
//! \details This is the public interface to most of the methods in the tree to
//! be used by the program.
struct CellRef {
private:
    //! the group of siblings this cell belongs to
    TreeGroup* _group;

    //! The index of the child within the sibling group.
    int _index;

public:
    //! the group of siblings this cell belongs to
    TreeGroup* group() const;

    //! The index of the child within the sibling group.
    int index() const;

    //! Returns the raw TreeCell data structure.
    TreeCell* treecell() const;

    //! Steps on to the next child in the current group
    void stepChild() {
        _index++;
    }

    //! Ctor for a cell reference.
    //! \details Cell references should be generated from a call to
    //!         \ref Tree::find(), or by traversal of the tree.
    //! \remark With debug enabled, this overload cannot be used to create the
    //!         special cases of the CellRef structure.
    CellRef(TreeGroup* gr, int ind);

    //! Ctor for an invalid cell ref
    CellRef();

    //! Ctor for a cell ref with flags
    //! \param flag The value put into the CellRef index. The group will be null.
    //!             A value of -1 will give the root cell. Any other value will
    //!             be the nth child of the root cell.
    CellRef(int flag);

    //! Checks if the reference is valid.
    bool isValid() const;

    //! Checks if the reference is the root cell special case.
    bool isRoot() const;

    //! Checks if the reference is to a boundary
    bool isBoundary() const;

    //! gets the boundary the cell represents
    int boundary() const;

    //! Reference to the parent cell
    CellRef parent() const;

    //! Returns the reference to the nth index child.
    //! \param n The child index. \copydetails CHILDREN_INDICES
    CellRef child(int n) const;

    //! Returns the nth child adjacent to a face.
    //! \param fc The face index (neighbour indexing).
    //!           \copydetails NEIGHBOUR_INDICES
    //! \param n The child cell reference. The order shouldn't be relied on.
    CellRef childOnFace(int fc, int n) const;

    //! Returns the child index adjacent to a face
    //! \param fc The face index (neighbour indexing).
    //!           \copydetails NEIGHBOUR_INDICES
    //! \param n The child index. The order shouldn't be relied on.
    static int childIndexOnFace(int fc, int n);

    //! Checks if the cell has children
    bool hasChildren() const;

    //! Checks if the cell has grand-children
    bool hasGrandChildren() const;

    //! Returns a pointer to the ortho containing the child cells
    TreeGroup* children() const;

    //! Gets a neighbouring cell
    //! \param n The neighbour index. \copydetails NEIGHBOUR_INDICES
    CellRef neighbour(int n) const;

    //! Gets a diagonal neighbouring cell
    //! \param n The diagonal index. These are indexed the same as vertices.
    //!          \copydetails VERTEX_INDICES
    CellRef diagonal(int n) const;

    //! Returns a convenient structure allowing iteration over neighbours
    tagNeighbours neighbours() const;

    //! Returns a convenient structure allowing iteration over poisson neighbours
    tagPoissonNeighbours poissonNeighbourhood() const;

    //! The orientation of the space filling curve at this cell.
    int orientation() const;

    //! The next cell in the space filling curve
    CellRef next() const;

    //! The previous cell in the space filling curve
    CellRef prev() const;

    //! Checks if the cell is the last child of this ortho in the space filling
    //! curve.
    bool isLastInGroup() const;

    //! Checks if the cell is the first child of this ortho in the space filling
    //! curve.
    bool isFirstInGroup() const;

    //! The processor this cell will be evaluated on
    node_t& rank();

    //! \copydoc CellRef::rank()
    node_t const& rank() const;

    //! The identifier of the cell
    //! \see TreeId
    ident_t id() const;

    //! The level of the cell in the tree.
    //! \remarks Zero is first level of children in the problem. Calling this
    //! method on the root cell special case will return -1.
    int level() const;

    //! The data associated with the cell.
    data_t& data();

    //! \copydoc CellRef::data()
    data_t const& data() const;

    //! Average values of all children (recursive)
    double avrChildren(fnData dt) const;

    //! Interpolate a "neighbours" data
    double ngbVal(int nb, fnData dt, double* ifBoundary = nullptr) const;

    // //! The data associated with a cell face.
    // facedata_t& facedata(int dir);
    // facedata_t const& facedata(int dir) const;

    //! The cell faces
    #ifdef TFTT_FACES
    FaceRef face(int dir) const;
    #endif

    #ifdef TFTT_VERTICES
    VertexRef vertex(int v) const;
    #endif

    //! The coordinate of the origin of the cell
    //! \details The origin is the lowest value.
    //! \param d The dimension to get.
    //! \remarks Use \ref CellRef::origins() instead if all dimensions are needed.
    double origin(int d) const;

    //! The coordinate of the origin of the cell
    //! \param[out] ret The origin vector.
    void origins(double ret[DIM]) const;

    //! The coordinate of the centre of the cell
    //! \param d The dimension to get.
    double centre(int d) const;

    //! The size of the cell
    //! \param d The dimension to get.
    //! \remarks Use \ref CellRef::sizes() instead if all dimensions are needed.
    double size(int d) const;

    //! The size of the cell (all dimensions)
    //! \param[out] ret The size vector.
    void sizes(double ret[DIM]) const;

    //! The location of a cell vertex
    //! \param v The vertex index. \copydetails VERTEX_INDICES
    //! \param d The dimension to get.
    double vertexPoint(int v, int d) const;

    //! The location of all the cell vertices
    //! \param[out] ret An array of vertex vectors.
    void vertices(double ret[1<<DIM][DIM]) const;

    //! Checks if the cell contains a point.
    //! \param pt A point in space.
    bool containsPoint(double pt[DIM]) const;

    bool operator==(const CellRef& rhs) const;
    bool operator!=(const CellRef& rhs) const;
    bool operator<(const CellRef& rhs) const;

    //! Access the data stored with the cell
    data_t* operator->();

    //! \copydoc CellRef::operator->()
    data_t const* operator->() const;

    //! A human readable cell identifier
    std::string path() const;

    //! Compares two cells
    //! \details The comparison is based on level, group then index in that order.
    //! \ref CellRef::parless should be used preferably unless level based
    //! comparison is necessary.
    struct less {
        bool operator()(const CellRef& a, const CellRef& b) {
            if (a.level() < b.level())
                return true;
            else if (a.level() > b.level())
                return false;
            else if (a.group() < b.group())
                return true;
            else if (a.group() > b.group())
                return false;
            else if (a.index() < b.index())
                return true;
            else
                return false;
        }
    };

    //! Compares two cells
    //! \details The comparison is based on the unique identifier of the cells.
    struct parless {
        bool operator()(const CellRef& a, const CellRef& b) {
            return a.id().id < b.id().id;
        }
    };
};


} // namespace tftt


//! Prints a convenient human readable description of the cell.
std::ostream& operator<<(std::ostream& os, const tftt::CellRef& cr);


#endif
