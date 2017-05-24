
#ifndef TFTT_PARALLEL_H
#define TFTT_PARALLEL_H


#include <boost/mpi/communicator.hpp>

#include "config.h"
#include "cellref.h"


namespace tftt {


//! Raw data sent over MPI to represent a cell
struct CellMin {
    ident_t id;
    node_t rank;
    data_t data;
};

//! Raw MPI data to signal a ghost rank change
struct GhostRankChange {
    ident_t id;
    node_t newRank;
};


//! Distribute cells across processors
//! \param n The number of nodes to distribute across
void distribute(int n);


//! Synchronises ghost data across the world
//! \param world The communicator object for the set of nodes
//! \remark Further documentation of the method can be found in the final report.
void syncGhosts(boost::mpi::communicator world);

//! Moves cells to the left/right nodes
//! \param world The communicator object for the set of nodes
//! \param left The number of cells to move to the left. If cells are being
//!             received this should be negative.
//! \param right The number of cells to move to the right. If cells are being
//!             received this should be negative.
//! \remark Further documentation of the method can be found in the final report.
void moveCells(boost::mpi::communicator world, int left, int right);


void adaptParSwBegin();
void adaptParSwPropogateVector(CellRef cl, uint8_t ref[2*DIM], uint8_t hold[2*DIM]);
void adaptParSwPropogateLevel(CellRef cl, int dir, int lvl, int offset=0);
void adaptParSwSetRefine(CellRef cl);


} // namespace tftt


#endif
