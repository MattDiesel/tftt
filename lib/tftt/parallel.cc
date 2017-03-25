
#include <iostream>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/request.hpp>

#include "tftt.h"

#include "parallel.h"

#include <unistd.h> // TODO: Remove


namespace mpi = boost::mpi;


namespace tftt {


void syncGhosts(mpi::communicator world)
{
    int wr = world.rank();
    int ws = world.size();

    int sz;

    mpi::request snt[ws];
    mpi::request rcv[ws];

    // Send all borders data
    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        sz = gtree.rawBorders[r].size();

        for (int b = 0; b < sz; b++) {
            gtree.borderData[r][b] = gtree.rawBorders[r][b].data();
        }

        snt[r] = world.isend(r, 0x100 | (wr*16) | r,
                             reinterpret_cast<int*>(gtree.borderData[r]), sz*sizeof(data_t)/sizeof(int));

        // std::cout << "[" << wr << "] Sent #" << (0x100 | (wr*16) | r) << " "
        //           << sz << " cells to " << r << "\n";
    }

    // Recieve all ghosts data
    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        sz = gtree.rawGhosts[r].size();

        rcv[r] = world.irecv(r, 0x100 | (r*16) | wr,
                             reinterpret_cast<int*>(gtree.ghostData[r]), sz*sizeof(data_t)/sizeof(int));

        // std::cout << "[" << wr << "] Recv #" << (0x100 | (r*16) | wr) << " "
        //           << sz << " cells from " << r << '\n';
    }


    // Wait for send/recv
    int mask = (1 << ws)-1;
    mask &= ~(1 << wr);
    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        if (!snt[r].test()) {
            // std::cout << "Wait for " << wr << "-" << r << " to send\n";
            snt[r].wait();
        }
    }

    mask = (1 << ws)-1;
    mask &= ~(1 << wr);
    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        if (!rcv[r].test()) {
            // std::cout << "Wait for " << r << "-" << wr << " to be recvd\n";
            rcv[r].wait();
        }
    }

    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        sz = gtree.rawGhosts[r].size();

        for (int g = 0; g < sz; g++) {
            gtree.rawGhosts[r][g].data() = gtree.ghostData[r][g];
        }
    }

    // std::cout << "Done.\n";
    world.barrier();
}


} // namespace tftt
