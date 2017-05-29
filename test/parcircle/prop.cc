
#include <iostream>
#include <fstream>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/operations.hpp>

#include "util/formatstring.h"
#include "util/pars.h"

#include "tftt/tftt.h"
#include "tftt/structure/treecell.h"

using namespace util; // for formatString
namespace mpi = boost::mpi;


int ITER = 0;

struct circle {

    double pos[2];
    double r;

    bool contains(double x, double y) {
        double dx = (x-pos[0]);
        double dy = (y-pos[1]);

        return dy*dy + dx*dx < r*r;
    }

    bool intersects(tftt::cell_t cl) {
        // Intersects if 1 vertex lies inside and 1 outside.

        bool in = false;
        bool out = false;
        for (int v = 0; v < 1<<DIM; v++) {
            if (contains(cl.vertexPoint(v, 0), cl.vertexPoint(v, 1))) {
                if (out) return true;
                in = true;
            }
            else {
                if (in) return true;
                out = true;
            }
        }

        return false;
    }
};

void printTwo2one(std::string fname, tftt::cell_t c)
{
    std::ofstream ofs(fname);

    tftt::cell_t nb;

    for (int dir = 0; dir < 4; dir++) {
        nb = c;
        for (int d = c.level(); d > 1; d--) {
            for (int p = 0; p < tftt::options.two2oneFlag; p++) {
                nb = nb.neighbour(dir);

                if (nb.isBoundary()) goto edge;

                tftt::drawCell(ofs, nb);
            }

            if (nb.level() > d)
                nb = nb.parent();
        }
    edge:
        nb=nb;
    }
}

bool areNeighbours(tftt::cell_t a, tftt::cell_t b)
{
    for (int nb = 0; nb < 2*DIM; nb++) {
        if (a.neighbour(nb) == b || b.neighbour(nb) == a) return true;
    }
    return false;
}

void miscChecks()
{
    std::set<uint64_t> seen;
    std::set<uint64_t> seenCurve;

    for (auto cl : tftt::leaves) {
        if (cl.isBoundary())
            throw std::runtime_error("Boundary in leaves");
        if (cl.hasChildren())
            throw std::runtime_error("Parent in leaves");

        auto result = seen.insert(cl.id().id);
        if (!result.second)
            throw std::runtime_error("Duplicate in leaf iterator");
    }

    tftt::cell_t prev;
    for (auto cl : tftt::curve) {
        if (cl.isBoundary())
            throw std::runtime_error("Boundary in curve");
        if (cl.hasChildren())
            throw std::runtime_error("Parent in curve");

        auto search = seen.find(cl.id().id);
        if (search == seen.end())
            throw std::runtime_error("In curve, not in leaves");
        seen.erase(cl.id().id);

        auto result = seenCurve.insert(cl.id().id);
        if (!result.second)
            throw std::runtime_error("Duplicate in curve iterator");

        if (prev.isValid()) {
            if (!areNeighbours(cl, prev))
                throw std::runtime_error("Consecutive cells in curve aren't neighbours");
            prev = cl;
        }
    }
    if (!seen.empty()) {
        throw std::runtime_error("Curve count != leaf count");
    }

    tftt::cell_t inside;
    for (int b = 0; b < 2*DIM; b++) {
        for (auto bcl : tftt::boundaryCells(b)) {
            if (bcl.hasChildren())
                throw std::runtime_error("Boundary leaf has children");
            if (!bcl.isBoundary())
                throw std::runtime_error("Boundary leaf isn't boundary");
            if (bcl.boundary() != b)
                throw std::runtime_error("Boundary leaf is on the wrong one");

            inside = bcl.neighbour(b ^ 1);
            if (inside.isBoundary())
                throw std::runtime_error("Boundary inside cell is still boundary");
            if (inside.level() != bcl.level())
                throw std::runtime_error("Boundary inside cell level mismatch");
            if (inside.hasChildren())
                throw std::runtime_error("Boundary inside cell has children");
        }
    }
}


int getAdaptVector(tftt::cell_t cl, int dir)
{
    tftt::TreeCell& tc = cl.group->cells[cl.index];
    return tc.adaptVector[dir];
}

uint32_t getAdaptVectors(tftt::cell_t cl)
{
    tftt::TreeCell& tc = cl.group->cells[cl.index];
    return *((uint32_t*)tc.adaptVector);
}

void setAdaptVectors(tftt::cell_t cl, uint32_t amount)
{
    tftt::TreeCell& tc = cl.group->cells[cl.index];
    for (int dir = 0; dir < 2*DIM; dir++) {
        tc.adaptVector[dir] = amount >> (dir*8);
        if (amount >> (dir*8)) std::cout << cl << " in direction " << (int)dir << " is of interest.\n";
    }
}


int getAdaptHoldVector(tftt::cell_t cl, int dir)
{
    tftt::TreeCell& tc = cl.group->cells[cl.index];
    return tc.adaptHoldVector[dir];
}

void project(float x, float y, float d, int dir, float& xo, float& yo)
{
    xo=x;
    yo=y;
    switch (dir) {
        case 0:
            xo -= d;
            break;
        case 1:
            xo += d;
            break;
        case 2:
            yo -= d;
            break;
        case 3:
            yo += d;
            break;
    }
}

void plotVect(std::ostream& os, float x, float y, int dir, float d)
{
    float x2,y2;
    os << x << " " << y << "\n";
    project(x, y, d, dir, x2, y2);
    os << x2 << " " << y2 << "\n\n";
}

void plotAdaptVectors(std::string fname)
{
    std::ofstream ofs(fname);

    int ad;
    for (int n = 0; n < 2; n++) {
        for (auto gh : tftt::gtree.rawBorders[n]) {
            for (int d = 0; d < 4; d++) {
                ad = getAdaptVector(gh, d);
                if (ad) {
                    plotVect(ofs, gh.centre(0), gh.centre(1), d, ad*gh.size(d >> 1));
                }
            }
        }
    }
}

void plotAdaptHoldVectors(std::string fname)
{
    std::ofstream ofs(fname);

    int ad;
    for (auto ghb : tftt::gtree.rawBorders) {
        for (auto gh : ghb) {
            for (int d = 0; d < 4; d++) {
                ad = getAdaptHoldVector(gh, d);
                if (ad) {
                    plotVect(ofs, gh.centre(0), gh.centre(1), d, ad*gh.size(d >> 1));
                }
            }
        }
    }
}

namespace tftt {

void syncVectors(mpi::communicator world)
{
    int wr = world.rank();
    int ws = world.size();

    int sz;

    mpi::request snt[ws];
    mpi::request rcv[ws];

    // Send all ghosts data
    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        sz = gtree.rawGhosts[r].size();

        for (int b = 0; b < sz; b++) {
            gtree.ghostAdaptVectors[r][b] = getAdaptVectors(gtree.rawGhosts[r][b]);
        }

        snt[r] = world.isend(r, 0x200 | (wr*16) | r,
                             reinterpret_cast<int*>(gtree.ghostAdaptVectors[r]), sz*sizeof(uint32_t)/sizeof(int));

        std::cout << "[" << wr << "] Sent #" << (0x200 | (wr*16) | r) << " "
                  << sz << " vectors to " << r << "\n";
    }

    // Recieve all borders data
    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        sz = gtree.rawBorders[r].size();

        rcv[r] = world.irecv(r, 0x200 | (r*16) | wr,
                             reinterpret_cast<int*>(gtree.borderAdaptVectors[r]), sz*sizeof(uint32_t)/sizeof(int));

        std::cout << "[" << wr << "] Recv #" << (0x200 | (r*16) | wr) << " "
                  << sz << " vectors from " << r << '\n';
    }


    // Wait for send/recv
    int mask = (1 << ws)-1;
    mask &= ~(1 << wr);
    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        if (!snt[r].test()) {
            std::cout << "Wait for " << wr << "-" << r << " to send\n";
            snt[r].wait();
        }
    }

    mask = (1 << ws)-1;
    mask &= ~(1 << wr);
    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        if (!rcv[r].test()) {
            std::cout << "Wait for " << r << "-" << wr << " to be recvd\n";
            rcv[r].wait();
        }
    }

    for (int r = 0; r < ws; r++) {
        if (r == wr) continue;

        sz = gtree.rawBorders[r].size();

        for (int g = 0; g < sz; g++) {
            setAdaptVectors(gtree.rawBorders[r][g], gtree.borderAdaptVectors[r][g]);
        }
    }

    // std::cout << "Done.\n";
    world.barrier();
}

} // namespace tftt


int main(int argc, char* argv[])
{
    mpi::environment env(argc, argv);
    mpi::communicator world;

    std::cout << "MPI rank=" << world.rank() << "/" << world.size() << '\n';


    // Defaults:
    int minDepth = 3;
    int maxDepth = 5;
    int iterations = 100;

    circle c;

    double startPos[2] {0.3, 0.3};
    double endPos[2] {0.7, 0.7};

    c.r = 0.2;
    c.pos[0] = startPos[0];
    c.pos[1] = startPos[1];

    tftt::options.two2oneFlag = 2;

    // Read from file:
    if (argc > 1) {
        try {
            getpars(argc, argv);

            tfetch("minDepth", minDepth);
            tfetch("maxDepth", maxDepth);
            tfetch("iterations", iterations);
            tfetch("circle.start[0]", startPos[0]);
            tfetch("circle.start[1]", startPos[1]);
            tfetch("circle.end[0]", endPos[0]);
            tfetch("circle.end[1]", endPos[1]);
            tfetch("circle.radius", c.r);
            tfetch("tftt.two2one", tftt::options.two2oneFlag);
        }
        catch (std::exception& e) {
            std::cout << "Error reading parameter file: " << e.what() << std::endl;
        }
    }
    else {
        if (!world.rank())
            std::cout << "Using default parameters" << std::endl;
    }

    if (!world.rank()) {

        std::cout << "Using Parameters: \n"
                  << "\tminDepth = " << minDepth << "\n"
                  << "\tmaxDepth = " << maxDepth << "\n"
                  << "\titerations = " << iterations << "\n"
                  << "\tcircle.start[0] = " << startPos[0] << "\n"
                  << "\tcircle.start[1] = " << startPos[1] << "\n"
                  << "\tcircle.end[0] = " << endPos[0] << "\n"
                  << "\tcircle.end[1] = " << endPos[1] << "\n"
                  << "\tcircle.radius = " << c.r << "\n"
                  << "\ttftt.two2one = " << tftt::options.two2oneFlag << std::endl;

        // Init tree to min depth
        tftt::init(1.0, 1.0);
        miscChecks();
        for (int d = 0; d < minDepth; d++) {
            for (auto& cl : tftt::leaves) {
                tftt::refine(cl);
            }
        }
        miscChecks();

        // Refine to circle.
        for (int d = minDepth; d < maxDepth; d++) {
            tftt::adaptSwBegin();
            for (auto& cl : tftt::leaves) {
                if (c.intersects(cl))
                    tftt::adaptSwSetRefine(cl);
            }
            tftt::adaptSwCommit();
        }
        miscChecks();
        for (auto& cl : tftt::leaves) {
            tftt::calcFaceCoefs(cl);
        }

        tftt::plot::mesh("parcircle/mesh.init.dat");

        tftt::distribute(world.size());
        tftt::saveParTree("parcircle.r{0}.tr", world.size());

        tftt::reset();
    }

    world.barrier();

    tftt::loadParTree(formatString("parcircle.r{0}.tr", world.rank()));

    tftt::plot::partialMesh(formatString("parcircle/mesh.r{0}.init.dat", world.rank()));
    for (int b = 0; b < world.size(); b++) {
        if (b == world.rank()) continue;

        tftt::plot::ghostMesh(formatString("parcircle/ghosts.r{0}.b{1}.init.dat", world.rank(), b), b);
        tftt::plot::borderMesh(formatString("parcircle/border.r{0}.b{1}.init.dat", world.rank(), b), b);

        std::cout << "[" << world.rank() << "] Ghosts to " << b
                  << ": " << tftt::gtree.ghosts[b].size() << "\n";
        std::cout << "[" << world.rank() << "] Borders to " << b
                  << ": " << tftt::gtree.borders[b].size() << "\n";
    }

    // What happens if we refine the cell at p
    tftt::adaptParSwBegin();

    for (auto& cl : tftt::leaves) {
        if (c.intersects(cl)) {
            if (cl.level() < maxDepth+1)
                tftt::adaptParSwSetRefine(cl);
            else if (cl.level() == maxDepth+1)
                tftt::adaptSwSetHoldRefined(cl);
        }
    }
    // for (auto& cl : tftt::leaforthos) {
    //     if (cl.level() >= minDepth) {
    //         tftt::adaptSwSetCoarsen(cl);
    //     }
    // }

    tftt::adaptSwCommit();

    tftt::syncVectors(world);


    plotAdaptVectors(formatString("parcircle/advec.r{0}.dat", world.rank()));
    plotAdaptHoldVectors(formatString("parcircle/hvec.r{0}.dat", world.rank()));


    tftt::adaptParSwBegin();
    uint8_t hold[4] = {0,0,0,0};
    for (int r = 0; r < world.size(); r++) {
        if (r == world.rank()) continue;

        for (auto bc : tftt::gtree.rawBorders[r]) {
            tftt::adaptParSwPropogateVector(bc, bc.group->cells[bc.index].adaptVector, hold);
        }
    }
    tftt::adaptSwCommit();

    tftt::plot::partialMesh(formatString("parcircle/mesh.r{0}.final.dat", world.rank()));
    for (int b = 0; b < world.size(); b++) {
        tftt::plot::ghostMesh(formatString("parcircle/ghosts.r{0}.b{1}.final.dat", world.rank(), b), b);
        tftt::plot::borderMesh(formatString("parcircle/border.r{0}.b{1}.final.dat", world.rank(), b), b);
    }

    return 0;
}
