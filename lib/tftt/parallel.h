
#ifndef TFTT_PARALLEL_H
#define TFTT_PARALLEL_H


#include <boost/mpi/communicator.hpp>

#include "tftt.h"
#include "tree.h"


// namespace boost {
// namespace mpi {

// template<> struct is_mpi_datatype<rt_data>
//     : public mpl::true_ { };

// }
// } // namespace boost::mpi


namespace tftt {


void syncGhosts(boost::mpi::communicator world);


} // namespace tftt


#endif
