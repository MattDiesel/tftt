
#ifndef TFTT_PUBLIC_H
#define TFTT_PUBLIC_H


namespace tftt {


// Constants for children, vertices, neighbours

//! Neighbours. P = positive dir, N = negative
//! Neighbours indices are made up of the direction (1 for positive) in the
//! lsb, and the dimension number in higher bits.
enum NEIGHBOUR_INDICES {
    NB_XN = 0,
    NB_XP,
    NB_YN,
    NB_YP,
    NB_ZN,
    NB_ZP
};


//! Vertices (2d)
//! Vertices are numbered in the same way children are numbered as seen by the
//! centre of the cell - using one bit per dimension, with a 1 in that position
//! representing the positive direction. Dimensions are ordered with the lowest
//! ranked dimension in the lsb, so {x,y,z} coordinates would be stored ZYX.
enum VERTEX_INDICES {
    VT_BL = 0,
    VT_BR = 1,
    VT_TL = 2,
    VT_TR = 3
};


//! Children (2d)
//! Children indexing is uses one bit per dimension, with a 1 in that position
//! representing the positive direction. Dimensions are ordered with the lowest
//! ranked dimension in the lsb, so {x,y,z} coordinates would be stored ZYX.
enum CHILDREN_INDICES {
    CH_BL = 0,
    CH_BR = 1,
    CH_TL = 2,
    CH_TR = 3
};


} // namespace tftt


#include "config.h"
#include "treeid.h"

#include "cellref.h"
#include "faceref.h"
#include "vertexref.h"


// bring utilities into tftt namespace
#include "iter/all.h"
#include "io/tfttio.h"
#include "io/serialize.h"
#include "tfttops.h"
#include "adapt.h"
#include "fttcore.h"

#include "structure/tree.h"

#ifndef TFTT_NOMPI
    #include "parallel.h"
#endif


namespace tftt {

// Miscellaneous functions that are mostly deprecated

//! Checks all the nearby cells satisfy the predicate
bool checkAround(cell_t cl, int dist, fnCheckCell check);

//! Finds any nearby cell that satisfies the predicate
bool findAround(cell_t cl, int dist, fnCheckCell check);

//! Interpolate a value on a cell face
//! \deprecated Cells should instead set the poisson coefficients and cell
//! references to acheive the effect.
double interpFace(cell_t cl, int fc, fnData dt);

//! Interpolates the value a child cell would have, if it existed.
//! \deprecated Cells should instead set the poisson coefficients and cell
//! references to acheive the effect.
double interpChild(cell_t cl, int ch, int forNb, fnData dt);

//! Interpolate the value at a vertex
//! \deprecated Vertex storage is currently being rewritten to store data on the
//! ALE grid, making this function redundant.
double interpALEVertex(cell_t cl, int v, fnData dt);


} // namespace tftt


#endif
