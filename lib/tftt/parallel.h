
#ifndef TFTT_PARALLEL_H
#define TFTT_PARALLEL_H


#include <boost/mpi/communicator.hpp>

#include "config.h"
#include "cellref.h"


// namespace boost {
// namespace mpi {

// template<> struct is_mpi_datatype<rt_data>
//     : public mpl::true_ { };

// }
// } // namespace boost::mpi


namespace tftt {


struct CellMin {
    ident_t id;
    node_t rank;
    data_t data;
};

struct GhostRankChange {
    ident_t id;
    node_t newRank;
};


//! Distribute cells across processors
void distribute(int n);

void syncGhosts(boost::mpi::communicator world);


void adaptParSwBegin();
void adaptParSwPropogateVector(CellRef cl, uint8_t ref[2*DIM], uint8_t hold[2*DIM]);
void adaptParSwPropogateLevel(CellRef cl, int dir, int lvl, int offset=0);
void adaptParSwSetRefine(CellRef cl);


void moveCells(boost::mpi::communicator world, int left, int right);


} // namespace tftt


#endif
