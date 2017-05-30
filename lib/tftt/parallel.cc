
#include <iostream>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/request.hpp>

#include "config.h"
#include "cellref.h"
#include "structure/tree.h"
#include "iter/all.h"
#include "adapt.h"
#include "fttcore.h"
#include "tfttops.h"

#include "parallel.h"


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


#ifdef TFTT_PARALLEL_REFPROP


void adaptParSwBegin()
{
    TreeCell* tc;
    for (auto cl : leaves) {
        tc = cl.treecell();
        tc->adaptFlags = AF_NoAction;

        tc = cl.parent().treecell();
        tc->adaptFlags = AF_NoAction;
    }

    for (auto ghb : gtree.rawGhosts) {
        for (auto gh : ghb) {

            TreeCell* tc = gh.treecell();
            for (int d = 0; d < 2*DIM; d++) {
                tc->adaptVector[d] = 0;
            }
        }
    }
}

void setAdaptVector(cell_t cl, int dir, int amount)
{
    TreeCell* tc = cl.treecell();

    if (amount > tc->adaptVector[dir]) {
        tc->adaptVector[dir] = amount;
    }
}

void setAdaptHoldVector(cell_t cl, int dir, int amount)
{
    TreeCell* tc = cl.treecell();

    if (amount > tc->adaptHoldVector[dir]) {
        tc->adaptHoldVector[dir] = amount;
    }
}

void adaptParSwPropogateLevel(cell_t cl, int dir, int lvl, int offset)
{
    int p = 0;

    if (offset) {
        p = options.two2oneFlag-offset;
    }

    for (int d = lvl; d > 1; d--) {
        for (; p < options.two2oneFlag; p++) {
            cl = cl.neighbour(dir);
            if (cl.isBoundary()) return;

            if (cl.hasChildren()) {
                // return; // I'm sure this should work.
            }
            else if (cl.level() == d) {
                adaptSwSetFlags(cl.parent(), AF_HoldRefined);
                adaptParSwPropogateLevel(cl, dir ^ 2, d-1);
                adaptParSwPropogateLevel(cl, dir ^ 3, d-1);
            }
            else if (cl.level() < d) {
                p++;
                adaptSwSetFlags(cl, AF_Refine);
                adaptParSwPropogateLevel(cl, dir ^ 2, d-1);
                adaptParSwPropogateLevel(cl, dir ^ 3, d-1);
            }

            if (cl.rank()) {
                if (cl.level() < d) {
                    setAdaptVector(cl, dir, options.two2oneFlag-p+1);
                    // std::cout << "R " << cl << "in " << dir << " = " << (options.two2oneFlag-p+1) << "\n";
                }
                else if (cl.level() == d) {
                    setAdaptHoldVector(cl, dir, options.two2oneFlag-p+1);
                    // std::cout << "h " << cl << "in " << dir << " = " << (options.two2oneFlag-p+1) << "\n";
                }
                // return; // Could be bad
            }
        }
        if (cl.level() == d) cl = cl.parent();
        p = 0;
    }
}

void adaptParSwPropogateVector(cell_t cl, uint8_t ref[2*DIM], uint8_t hold[2*DIM])
{
    for (int d = 0; d < 2*DIM; d++) {
        if (ref[d]) {
            std::cout << "In dir " << d << " = " << (int)ref[d] << "\n";
            adaptSwSetFlags(cl, AF_Refine);
            adaptParSwPropogateLevel(cl, d, cl.level(), ref[d]-1);
        }
        else if (hold[d]) {
            adaptParSwPropogateLevel(cl, d, cl.level(), hold[d]-1);
        }
    }
}

void adaptParSwSetRefine(cell_t cl)
{
    if (cl.hasChildren()) {
        throw std::logic_error("Cell already refined.");
    }

    adaptSwSetFlags(cl, AF_Refine);
    for (int d = 0; d < 4; d++) {
        adaptParSwPropogateLevel(cl, d, cl.level());
    }
}


#endif


void addNeighb(std::set<cell_t>& ghosts, cell_t cl, node_t node)
{
    // Ghosts are any cells required by poisson coefficients

    TreeCell* tc = cl.treecell();

    for (int n = 0; n < tc->poisNgbC; n++) {
        if (tc->poisNgb[n].rank() != node && !tc->poisNgb[n].isBoundary()) {
            ghosts.insert(tc->poisNgb[n]);
        }
    }
}

void moveCells(boost::mpi::communicator world, int left, int right)
{
    // We know how many cells
    // We now need to know how many ghosts
    uint16_t fromLeft, fromRight;
    std::set<cell_t> toSendLeft;
    std::set<cell_t> toSendRight;
    CellMin* sendLeftVector = nullptr;
    CellMin* sendRightVector = nullptr;
    CellMin* recvLeftVector = nullptr;
    CellMin* recvRightVector = nullptr;

    std::vector<std::vector<GhostRankChange> > rankChanges;
    for (int r = 0; r < world.size(); r++) {
        rankChanges.push_back(decltype(rankChanges)::value_type());
    }

    // Update cell ranks
    // Must be done first, as any neighbours need to have the new up to date thread info
    if (left > 0) {
        tagActiveCurve::curve_iterator it = activecurve.begin();

        for (int i = 0; i < left; i++) {
            it->rank() = world.rank()-1;

            // Work out the effect of this rank change
            for (int b = 0; b < (int)gtree.borders.size(); b++) {
                if (b == world.rank() || b == world.rank()-1) continue;

                auto search = gtree.borders[b].find(*it);
                if (search != gtree.borders[b].end()) {
                    // std::cout << *it << " is a border cell to " << b << "\n";
                    rankChanges[b].push_back({it->id(), node_t(world.rank()-1)});
                }
            }

            it++;
        }
    }
    if (right > 0) {
        tagActiveCurve::curve_iterator it = activecurve.rbegin();

        for (int i = 0; i < right; i++) {
            it->rank() = world.rank()+1;

            // Work out the effect of this rank change
            for (int b = 0; b < (int)gtree.borders.size(); b++) {
                if (b == world.rank() || b == world.rank()+1) continue;

                auto search = gtree.borders[b].find(*it);
                if (search != gtree.borders[b].end()) {
                    // std::cout << *it << " is a border cell to " << b << "\n";
                    rankChanges[b].push_back({it->id(), node_t(world.rank()+1)});
                }
            }

            it--;
        }
    }


    // Notify rank changes
    for (int r = 0; r < world.size(); r++) {
        if (r == world.rank()) continue;

        world.isend(r, 0x500 | (world.rank() << 8) | r, (int)rankChanges[r].size());
    }
    std::vector<int> numChanges(world.size());
    for (int r = 0; r < world.size(); r++) {
        if (r == world.rank()) continue;

        world.recv(r, 0x500 | (r << 8) | world.rank(), numChanges[r]);
    }

    // Send changes
    for (int r = 0; r < world.size(); r++) {
        if (r == world.rank()) continue;
        if (rankChanges[r].size() == 0) continue;

        world.isend(r, 0x600 | (world.rank() << 8) | r,
                    reinterpret_cast<int*>(rankChanges[r].data()),
                    rankChanges[r].size()*sizeof(GhostRankChange)/sizeof(int));
    }
    // Recieve changes
    for (int r = 0; r < world.size(); r++) {
        if (r == world.rank()) continue;
        if (numChanges[r] == 0) continue;

        GhostRankChange* rc = new GhostRankChange[numChanges[r]];

        world.recv(r, 0x600 | (r << 8) | world.rank(),
                   reinterpret_cast<int*>(rc),
                   numChanges[r]*sizeof(*rc)/sizeof(int));

        for (int c = 0; c < numChanges[r]; c++) {
            cell_t cl = tftt::find(rc[c].id);
            cl.rank() = rc[c].newRank;

            // std::cout << "[" << world.rank() << "] Ghost " << cl << " rank changed to " << cl.rank() << "\n";
        }

        delete[] rc;
    }



    // Work out how many ghosts to send as well
    if (left > 0) {
        tagActiveCurve::curve_iterator it = activecurve.begin();

        for (int i = 0; i < left; i++) {
            toSendLeft.insert(*it);

            addNeighb(toSendLeft, *it, world.rank()-1);

            it++;
        }

        // std::cout << "[" << world.rank() << "] Sending " << left
        //           << " (" << toSendLeft.size() << ") left.\n";

        world.isend(world.rank()-1, 0x301 | (world.rank() << 8), (uint16_t)toSendLeft.size());
    }
    if (right > 0) {
        tagActiveCurve::curve_iterator it = activecurve.rbegin();

        for (int i = 0; i < right; i++) {
            toSendRight.insert(*it);

            addNeighb(toSendRight, *it, world.rank()+1);

            it--;
        }

        // std::cout << "[" << world.rank() << "] Sending " << right
        //           << " (" << toSendRight.size() << ") right.\n";

        world.isend(world.rank()+1, 0x302 | (world.rank() << 8), (uint16_t)toSendRight.size());
    }


    // Prepare send vectors
    if (left > 0) {
        sendLeftVector = new CellMin[toSendLeft.size()];
        int x = 0;
        for (auto cl : toSendLeft) {
            sendLeftVector[x].id = cl.id();
            sendLeftVector[x].rank = cl.rank();
            sendLeftVector[x].data = cl.data();

            x++;
        }
    }
    if (right > 0) {
        sendRightVector = new CellMin[toSendRight.size()];
        int x = 0;
        for (auto cl : toSendRight) {
            sendRightVector[x].id = cl.id();
            sendRightVector[x].rank = cl.rank();
            sendRightVector[x].data = cl.data();

            x++;
        }
    }


    // Recieve true recv amounts
    if (left < 0) {
        // Recv from left
        world.recv(world.rank()-1, 0x302 | ((world.rank()-1) << 8), fromLeft);

        // std::cout << "[" << world.rank() << "] Recv " << fromLeft << " from left\n";
    }
    if (right < 0) {
        // Recv from right
        world.recv(world.rank()+1, 0x301 | ((world.rank()+1) << 8), fromRight);

        // std::cout << "[" << world.rank() << "] Recv " << fromRight << " from right\n";
    }


    // Send data
    if (left > 0) {
        world.isend(world.rank()-1, 0x401 | (world.rank() << 8),
                    reinterpret_cast<int*>(sendLeftVector),
                    toSendLeft.size()*sizeof(*sendLeftVector)/sizeof(int));
    }
    if (right > 0) {
        world.isend(world.rank()+1, 0x402 | (world.rank() << 8),
                    reinterpret_cast<int*>(sendRightVector),
                    toSendRight.size()*sizeof(*sendRightVector)/sizeof(int));
    }


    // Recieve Data
    if (left < 0) {
        recvLeftVector = new CellMin[fromLeft];

        world.recv(world.rank()-1, 0x402 | ((world.rank()-1) << 8),
                   reinterpret_cast<int*>(recvLeftVector),
                   fromLeft*sizeof(CellMin)/sizeof(int));
    }
    if (right < 0) {
        recvRightVector = new CellMin[fromRight];

        world.recv(world.rank()+1, 0x401 | ((world.rank()+1) << 8),
                   reinterpret_cast<int*>(recvRightVector),
                   fromRight*sizeof(CellMin)/sizeof(int));
    }


    // Insert into ourselves
    cell_t newcl;
    if (left < 0) {
        for (int i = 0; i < fromLeft; i++) {
            newcl = tftt::insert(recvLeftVector[i].id);
            newcl.rank() = recvLeftVector[i].rank;
            newcl.data() = recvLeftVector[i].data;

            // std::cout << "[" << world.rank() << "] Added>> " << newcl
            //           << " rank " << (int)newcl.rank() << "\n";
        }
        for (int i = 0; i < -left; i++) {
            gtree.firstActive = gtree.firstActive.prev();
        }
    }
    if (right < 0) {
        for (int i = 0; i < fromRight; i++) {
            newcl = tftt::insert(recvRightVector[i].id);
            newcl.rank() = recvRightVector[i].rank;
            newcl.data() = recvRightVector[i].data;

            // std::cout << "[" << world.rank() << "] Added<< " << newcl
            //           << " rank " << (int)newcl.rank() << "\n";
        }
        for (int i = 0; i < -right; i++) {
            gtree.lastActive = gtree.lastActive.next();
        }
    }

    // Calculate poisson coefficients for new cells
    if (left < 0) {
        newcl = gtree.firstActive;
        for (int i = 0; i < -left; i++) {
            tftt::calcFaceCoefs(newcl);
            newcl = newcl.next();
        }
    }
    if (right < 0) {
        newcl = gtree.lastActive;
        for (int i = 0; i < -right; i++) {
            tftt::calcFaceCoefs(newcl);
            newcl = newcl.prev();
        }
    }

    // Roll back unused cells
    // These may be coarsened later.
    if (left > 0) {
        // Update start pointer
        for (int i = 0; i < left; i++) {
            gtree.firstActive = gtree.firstActive.next();
        }
    }
    if (right > 0) {
        for (int i = 0; i < right; i++) {
            gtree.lastActive = gtree.lastActive.prev();
        }
    }

    // May be required?
    for (auto cl : activecurve) {
        cl.rank() = world.rank();
    }

    // Update ghost cells
    for (int r = 0; r < world.size(); r++) {
        if (r == world.rank()) continue;
        gtree.ghosts[r].clear();
        gtree.borders[r].clear();
    }
    for (auto cl : activecurve) {
        TreeCell* tc = cl.treecell();

        for (int p = 0; p < tc->poisNgbC; p++) {
            if (tc->poisNgb[p].rank() != world.rank() && !tc->poisNgb[p].isBoundary()) {
                gtree.ghosts[tc->poisNgb[p].rank()].insert(tc->poisNgb[p]);
            }
        }
    }

    // Fill raw ghost
    for (int r = 0; r < (int)gtree.ghosts.size(); r++) {
        if (r == world.rank()) continue;

        gtree.rawGhosts[r].clear();

        for (auto& gh : gtree.ghosts[r]) {
            gtree.rawGhosts[r].push_back(gh);

            // Border cells are any required by the ghosts poisson coefs.
            calcFaceCoefs(gh);

            TreeCell* tc = gh.treecell();

            for (int p = 0; p < tc->poisNgbC; p++) {
                if (tc->poisNgb[p].rank() == gtree.rank) {
                    gtree.borders[gh.rank()].insert(tc->poisNgb[p]);
                }
            }
        }

        // std::cout << "[" << world.rank() << "] Now has " << gtree.ghosts[r].size()
        //           << " ghosts to " << r << "\n";
        // std::cout << "[" << world.rank() << "] Now has " << gtree.borders[r].size()
        //           << " borders to " << r << "\n";
    }

    // Fill raw
    for (int r = 0; r < (int)gtree.borders.size(); r++) {
        if (r == world.rank()) continue;
        gtree.rawBorders[r].clear();

        for (auto& br : gtree.borders[r]) {
            gtree.rawBorders[r].push_back(br);
        }
    }

    // Fill empty data array
    // Todo: Work out if there's been a change.
    size_t szg, szb;
    for (int r = 0; r < (int)gtree.ghosts.size(); r++) {
        if (r == world.rank()) continue;
        szg = gtree.ghosts[r].size();
        szb = gtree.borders[r].size();

        if (gtree.ghostData[r]) {
            delete[] gtree.ghostData[r];
            delete[] gtree.ghostAdaptVectors[r];
            gtree.ghostData[r] = nullptr;
            gtree.ghostAdaptVectors[r] = nullptr;
        }
        if (gtree.borderData[r]) {
            delete[] gtree.borderData[r];
            delete[] gtree.borderAdaptVectors[r];
            gtree.borderData[r] = nullptr;
            gtree.borderAdaptVectors[r] = nullptr;
        }

        if (szg) {
            gtree.ghostData[r] = new data_t[szg];
            gtree.ghostAdaptVectors[r] = new uint32_t[szg];
        }

        if (szb) {
            gtree.borderData[r] = new data_t[szb];
            gtree.borderAdaptVectors[r] = new uint32_t[szb];
        }
    }


    delete[] sendLeftVector;
    delete[] sendRightVector;
    delete[] recvLeftVector;
    delete[] recvRightVector;
}


} // namespace tftt
