
#include <iostream>
#include <fstream>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/operations.hpp>

#include "util/formatstring.h"
#include "util/pars.h"

#include "tftt/tftt.h"
#include "tftt/tree.h"

#include <vector>

namespace mpi = boost::mpi;
using namespace util; // for formatString


struct circle {

    double pos[2];
    double r;

    // @brief Checks if the circle contains a point.
    bool contains(double x, double y) {
        double dx = (x-pos[0]);
        double dy = (y-pos[1]);

        return dy*dy + dx*dx < r*r;
    }

    double cc(tftt::cell_t cl) {
        if (intersects(cl)) {
            // Naive method of computing VoF at this point
            double x1 = cl.origin(0), x2 = cl.origin(0) + cl.size(0);
            double y1 = cl.origin(1), y2 = cl.origin(1) + cl.size(1);
            double dx = cl.size(0) * 0.01, dy = cl.size(1)*0.01;

            double ret = 0.0;
            int count = 0;
            for (double x = x1; x <= x2; x += dx) {
                for (double y = y1; y <= y2; y += dy) {
                    if (contains(x, y))
                        ret += 1.0;
                    count++;
                }
            }

            ret /= count;
            return ret;
        }
        else {
            if (contains(cl.origin(0), cl.origin(1)))
                return 1.0;
            else
                return 0.0;
        }
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



int ITER = 0;
circle c;


double fn(tftt::cell_t& cl)
{
    return -cl->cc;
}

double dtPc(tftt::data_t& dt)
{
    return dt.P;
}

double dtPcn(tftt::data_t& dt)
{
    return -dt.P;
}

double& dtP(tftt::data_t& dt)
{
    return dt.P;
}


int main(int argc, char* argv[])
{
    mpi::environment env(argc, argv);
    mpi::communicator world;

    std::cout << "MPI rank=" << world.rank() << "/" << world.size() << '\n';

    // Defaults:
    int minDepth = 4;
    int maxDepth = 6;
    int iterations = 100;
    int plotEvery = 1;
    int printEvery = 1;
    int residEvery = 1;

    double initialValue = 0.0;

    double omega = 1.3;

    double startPos[2] {0.3, 0.3};
    double endPos[2] {0.7, 0.7};

    c.r = 0.2;
    c.pos[0] = startPos[0];
    c.pos[1] = startPos[1];

    tftt::options.two2oneFlag = 2;

    for (int b = 0; b < 2*DIM; b++) {
        tftt::gtree.isNeuman[b] = false;
        tftt::gtree.dirichletValue[b] = 0.0;
    }

    // Read from file:
    if (argc > 1) {
        try {
            getpars(argc, argv);

            tfetch("plotEvery", plotEvery);
            tfetch("printEvery", printEvery);
            tfetch("residEvery", residEvery);
            tfetch("minDepth", minDepth);
            tfetch("maxDepth", maxDepth);
            tfetch("omega", omega);
            tfetch("initialValue", initialValue);
            tfetch("iterations", iterations);
            tfetch("circle.start[0]", startPos[0]);
            tfetch("circle.start[1]", startPos[1]);
            tfetch("circle.end[0]", endPos[0]);
            tfetch("circle.end[1]", endPos[1]);
            tfetch("circle.radius", c.r);
            tfetch("tftt.two2one", tftt::options.two2oneFlag);

            tfetch("neuman[0]", tftt::gtree.isNeuman[0]);
            tfetch("dirichlet[0]", tftt::gtree.dirichletValue[0]);
            tfetch("neuman[1]", tftt::gtree.isNeuman[1]);
            tfetch("dirichlet[1]", tftt::gtree.dirichletValue[1]);
            tfetch("neuman[2]", tftt::gtree.isNeuman[2]);
            tfetch("dirichlet[2]", tftt::gtree.dirichletValue[2]);
            tfetch("neuman[3]", tftt::gtree.isNeuman[3]);
            tfetch("dirichlet[3]", tftt::gtree.dirichletValue[3]);
        }
        catch (std::exception& e) {
            std::cout << "Error reading parameter file: " << e.what() << std::endl;
        }
    }
    else {
        if (world.rank() == 0)
            std::cout << "Using default parameters" << std::endl;
    }
    if (world.rank() == 0) {

        std::cout << "Using Parameters: \n"
                  << "\tminDepth = " << minDepth << "\n"
                  << "\tmaxDepth = " << maxDepth << "\n"
                  << "\titerations = " << iterations << "\n"
                  << "\tresidEvery = " << residEvery << "\n"
                  << "\tcircle.start[0] = " << startPos[0] << "\n"
                  << "\tcircle.start[1] = " << startPos[1] << "\n"
                  << "\tcircle.end[0] = " << endPos[0] << "\n"
                  << "\tcircle.end[1] = " << endPos[1] << "\n"
                  << "\tcircle.radius = " << c.r << "\n"
                  << "\ttftt.two2one = " << tftt::options.two2oneFlag << std::endl;

        std::cout << "Boundary Conditions: \n";
        for (int b = 0; b < 4; b++) {
            std::cout << "\t[" << b << "] ";

            if (tftt::gtree.isNeuman[b])
                std::cout << "Neuman\n";
            else
                std::cout << "Dirichlet = " << tftt::gtree.dirichletValue[b] << "\n";
        }

        // Init tree to min depth
        tftt::init(1.0, 1.0);

        for (int d = 0; d < minDepth; d++) {
            for (auto& cl : tftt::leaves) {
                tftt::refine(cl);
            }
        }

        // Refine to circle.
        for (int d = minDepth; d < maxDepth; d++) {
            tftt::adaptBegin();

            for (auto& cl : tftt::leaves) {
                if (c.intersects(cl))
                    tftt::adaptAdd(cl);
            }

            tftt::adaptCommit();
        }

        tftt::drawMesh("parp/cmesh.init.dat");
        tftt::drawCurve("parp/chilb.init.dat");

        for (auto& cl : tftt::leaves) {
            tftt::calcFaceCoefs(cl);
        }

        tftt::distribute(world.size());
        tftt::splitToDisk("parp/r{0}.ltr");

        tftt::reset();
    }

    world.barrier();

    tftt::loadTree(formatString("parp/r{0}.ltr", world.rank()), world.rank());
    // tftt::loadTree(formatString("parp/r{0}.ltr", 0), 0);

    // Initial Conditions for cells
    for (auto& cl : tftt::activecurve) {
        cl->cc = c.cc(cl);
        cl->P = initialValue;
        tftt::calcFaceCoefs(cl);
    }

    for (int b = 0; b < 4; b++) {
        if (!tftt::gtree.isNeuman[b]) {
            for (auto& cl : tftt::boundaryCells(b)) {
                cl->P = tftt::gtree.dirichletValue[b];
            }
        }
    }

    tftt::drawMesh(formatString("parp/cmesh.r{0}.init.dat", world.rank()));
    tftt::drawCurve(formatString("parp/chilb.r{0}.init.dat", world.rank()));
    tftt::drawPartialMesh(formatString("parp/mesh.r{0}.init.dat", world.rank()));
    tftt::drawPartialCurve(formatString("parp/hilb.r{0}.init.dat", world.rank()));
    tftt::drawBoundaries(formatString("parp/bound.r{0}.init.dat", world.rank()));

    for (int n = 0; n < world.size(); n++) {
        tftt::drawGhosts(formatString("parp/ghosts.r{0}.b{1}.init.dat", world.rank(), n), n);
        tftt::drawBorder(formatString("parp/border.r{0}.b{1}.init.dat", world.rank(), n), n);
    }

    std::cout << "[" << world.rank() << "] Loaded tree of "
              << tftt::gtree.cactive << " cells (" << tftt::gtree.ccells << " total).\n";

    // world.barrier();
    tftt::syncGhosts(world);

    tftt::cell_t mx;

    double resid;
    double maxResid;
    for (ITER = 0; ITER < iterations; ITER++) {

        // mx = tftt::max(dtPc);
        // std::cerr << "Max P = " << mx->P << " @ " << mx << "\n";

        // mx = tftt::max(dtPcn);
        // std::cerr << "Min P = " << mx->P << " @ " << mx << "\n";

        tftt::relax(omega, dtP, fn);

        for (int b = 0; b < 4; b++) {
            if (!tftt::gtree.isNeuman[b]) {
                for (auto& cl : tftt::boundaryCells(b)) {
                    cl->P = tftt::gtree.dirichletValue[b];
                }
            }
        }

        tftt::syncGhosts(world);

        if (ITER % residEvery == 0) {
            resid = tftt::resid(dtP, fn);

            if (world.rank() == 0) {
                mpi::reduce(world, resid, maxResid, mpi::maximum<double>(), 0);
            }
            else {
                mpi::reduce(world, resid, mpi::maximum<double>(), 0);
            }
        }

        if (ITER % printEvery == 0 && world.rank() == 0) {
            std::cerr << "Iteration " << ITER << ", Resid: " << maxResid << '\n';
        }

        if (ITER % plotEvery == 0) {
            tftt::plotMatrix(formatString("parp/P.{0}.r{1}.dat", ITER, world.rank()), [](tftt::data_t& dt) {
                return dt.P;
            });

            tftt::plotMatrix(formatString("parp/res.{0}.r{1}.dat", ITER, world.rank()), [](tftt::data_t& dt) {
                return dt.res;
            });
        }
    }

    tftt::plotMatrix(formatString("parp/P.final.r{0}.dat", world.rank()), [](tftt::data_t& dt) {
        return dt.P;
    });
    tftt::plotMatrix(formatString("parp/res.final.r{0}.dat", world.rank()), [](tftt::data_t& dt) {
        return dt.res;
    });

    return 0;
}
