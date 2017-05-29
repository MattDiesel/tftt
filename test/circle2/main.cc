
#include <iostream>
#include <fstream>

#include "util/formatstring.h"
#include "util/pars.h"

#define TFTT_NOMPI

#include "tftt/tftt.h"
#include "tftt/tree.h"

using namespace util; // for formatString


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


int main(int argc, char* argv[])
{
    // Defaults:
    int minDepth = 2;
    int maxDepth = 6;
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
        std::cout << "Using default parameters" << std::endl;
    }

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

    tftt::plot::mesh("circle2data/mesh.min.dat");

    // Refine to circle.
    for (int d = minDepth; d < maxDepth; d++) {
        tftt::adaptSwBegin();

        for (auto& cl : tftt::leaves) {
            if (c.intersects(cl))
                tftt::adaptSwSetRefine(cl);
        }

        tftt::adaptSwCommit();

        tftt::plot::mesh(formatString("circle2data/mesh.init.{0}.dat", d-minDepth));
    }

    tftt::plot::mesh("circle2data/mesh.init.dat");

    tftt::cell_t tmp;
    for (ITER = 1; ITER <= iterations; ITER++) {
        c.pos[0] = startPos[0] + (endPos[0] - startPos[0])/iterations * ITER;
        c.pos[1] = startPos[1] + (endPos[1] - startPos[1])/iterations * ITER;

        std::cout << "Iteration " << ITER << "\n\tCircle at: " << c.pos[0] << "," << c.pos[1] << "\n";

        miscChecks();

        tftt::adaptSwBegin();
        for (auto& cl : tftt::leaves) {
            if (c.intersects(cl)) {
                if (cl.level() < maxDepth)
                    tftt::adaptSwSetRefine(cl);
                else if (cl.level() == maxDepth)
                    tftt::adaptSwSetHoldRefined(cl);
            }
        }
        for (auto& cl : tftt::leaforthos) {
            if (cl.level() >= minDepth) {
                tftt::adaptSwSetCoarsen(cl);
            }
        }
        tftt::adaptSwCommit();

        std::cout << "\tCell Count: " << tftt::gtree.ccells << "\n";

        tftt::plot::mesh(formatString("circle2data/mesh.{0}.dat", ITER));
    }

    return 0;
}
