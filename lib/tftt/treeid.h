
#ifndef LTREE_IDENT_H
#define LTREE_IDENT_H


#include <type_traits>
#include <cstdint>
#include <ostream>


namespace tftt {


//! A Cell Identifier in a TFT Tree
//! Encapsulates identifier encoding of path and tree level.
//! 
//! The storage format is variable sized, and can be arbitrarily large (provided a datatype exists).
//!
//! The most significant byte encodes the level, with the first group of cells being level zero. This is always fixed at
//! a single byte regardless of datatype size.
//!
//! The remaining bits encode the path, in tuples of bits depending on the dimension of the tree (1 bit per dimension),
//! each tuple is the index of the child, as defined in treecellorder.h. 
//!
//! This encoding scheme guarentees that an identifier for a cell will be unique within the tree. The identifier is also
//! sufficient to know the cells location within the tree, making it ideal for use across processing unit boundaries. 
template<
            typename T,         //!< The base type of the identifier. Must be an unsigned arithmetic type.
            unsigned int DIM    //!< The dimension of the tree. 
        >
struct TreeId {
    // Template Parameter Checks
    static_assert(std::is_unsigned<T>(), "ID Type must be an unsigned integral type.");
    static_assert(sizeof(T) >= 2, "ID Type must be more than one byte.");
    static_assert(DIM != 0, "Dimension must be non-zero.");


    // Constants used in bit operations

    static constexpr T LEVELSHIFT = ((sizeof(T)-1)*8);              //!< Index of the LSB of the level.
    static constexpr T LEVELMASK = (T(0xFF) << LEVELSHIFT);         //!< Mask of the level byte.
    static constexpr T LEVELINC = (T(0x01) << LEVELSHIFT);          //!< A level value of 1
    static constexpr T PATHMASK = (~(T(-1) << DIM));                //!< Mask of the bottom level path bits
    static constexpr T FULLPATHMASK = (~LEVELMASK);                 //!< Mask of all path bytes.

    //! The integral identifier value.
    T id;

    //! The level of the cell within the tree
    inline unsigned char level() const {
        return id >> LEVELSHIFT;
    }

    //! The encoded full path of the tree.
    inline T path() const {
        return id & FULLPATHMASK;
    }

    //! Returns the identifier of the parent cell
    inline TreeId parent() const {
        return (path() >> DIM) | ((id & LEVELMASK) - LEVELINC);
    }

    //! Returns the orthant (quadrant in 2d, octant in 3d) of the bottom level cell
    inline unsigned int orthant() const {
        return id & PATHMASK;
    }

    //! Returns the orthant of the nth level cell.
    inline unsigned int orthant(
                int n                   //!< The level to retrieve the orthant of, the top cell in the tree is level zero.
            ) const {
        return (id >> ((level() - n) * DIM)) & PATHMASK;
    }

    //! Returns the id of a child cell
    inline TreeId child(
                unsigned int orth       //!< The orthant of the child cell
            ) const {
        return firstchild().id | orth;
    }

    //! Returns the id of the zeroth child, also used as the id of a group.
    inline TreeId firstchild() const {
        return (path() << DIM) | ((id & LEVELMASK) + LEVELINC);
    }

    //! Test two identifiers for equality
    bool operator==(
                TreeId<T, DIM> rhs        //!< The second identifier to be compared.
            ) const {
        return id == rhs.id;
    }

    //! Ctor for the top level group
    //! Note: This will be the ID of the first child of the top level group.
    TreeId() : id(0) {}

    //! Ctor for an identifier with a given raw value.
    TreeId(T n) : id(n) {}
    // operator T() {
    //  return id;
    // }
};


// Type definitions of common identifier types.

// typedef TreeId<uint32_t, 2> qtid_t;           //!< 32 bit Quadtree identifier
// typedef TreeId<uint64_t, 2> qtidl_t;          //!< 64 bit Quadtree identifier
// typedef TreeId<uint32_t, 3> otid_t;           //!< 32 bit Octree identifier
// typedef TreeId<uint64_t, 3> otidl_t;          //!< 64 bit Octree identifier


} // namespace tftt


//! Writes the full path to a cell to a stream.
//! The path is human readable, the resulting string will be unique within the tree.
template<typename T, unsigned int DIM>
std::ostream& operator<<(
            std::ostream& os,               //!< The output stream to write to
            tftt::TreeId<T,DIM> id          //!< The identifier to write
        ) {
    for (int i = 0; i < id.level(); i++) {
        os << id.orthant(i) << "-";
    }
    os << id.orthant();
    return os;
}


#endif
