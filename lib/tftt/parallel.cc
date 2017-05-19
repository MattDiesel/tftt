
#include <iostream>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/request.hpp>

#include "config.h"
#include "tree.h"
#include "leaves.h"

#include "parallel.h"

#include <unistd.h> // TODO: Remove


namespace mpi = boost::mpi;


namespace tftt {


void distribute(int n)
{
    int cpernode = (gtree.ccells / n)+1;

    int cc = cpernode;
    int node = 0;
    for (auto& cl : curve) {
        if (!--cc) {
            cc = cpernode;
            node++;
        }

        cl.rank() = node;
    }
}


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


// void adaptParSwPropogateLevel(CellRef cl, int dir, int lvl)
// {
//     for (int d = lvl; d > 1; d--) {
//         for (int p = 0; p < options.two2oneFlag; p++) {
//             cl = cl.neighbour(dir);
//             if (cl.isBoundary()) return;

//             if (cl.hasChildren()) {
//                 // return; // I'm sure this should work.
//             }
//             else if (cl.level() == d) {
//                 adaptSwSetFlags(cl.parent(), AF_HoldRefined);
//                 adaptSwPropogateLevel(cl, dir ^ 2, d-1);
//                 adaptSwPropogateLevel(cl, dir ^ 3, d-1);
//             }
//             else if (cl.level() < d) {
//                 p++;
//                 adaptSwSetFlags(cl, AF_Refine);
//                 adaptSwPropogateLevel(cl, dir ^ 2, d-1);
//                 adaptSwPropogateLevel(cl, dir ^ 3, d-1);
//             }

//             if (cl.rank()) {
//                 if (cl.level() <= d)
//                     std::cout << "Propagate through ghost " << cl << " - " << dir << " " << d << "\n";
//             }
//         }
//         if (cl.level() == d) cl = cl.parent();
//     }
// }

// void adaptParSwSetRefine(CellRef cl)
// {
//     if (cl.hasChildren()) {
//         throw std::logic_error("Cell already refined.");
//     }

//     adaptSwSetFlags(cl, AF_Refine);
//     for (int d = 0; d < 4; d++) {
//         adaptParSwPropogateLevel(cl, d, cl.level());
//     }
// }


} // namespace tftt
