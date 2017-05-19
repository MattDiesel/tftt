
#ifndef TFTT_PARALLEL_H
#define TFTT_PARALLEL_H


#include <boost/mpi/communicator.hpp>

#include "config.h"


// namespace boost {
// namespace mpi {

// template<> struct is_mpi_datatype<rt_data>
//     : public mpl::true_ { };

// }
// } // namespace boost::mpi


namespace tftt {


//! Distribute cells across processors
void distribute(int n);

void syncGhosts(boost::mpi::communicator world);

// void adaptParSwPropogateLevel(CellRef cl, int dir, int lvl);
// void adaptParSwSetRefine(CellRef cl);


} // namespace tftt


#endif
