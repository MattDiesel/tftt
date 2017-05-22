
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



int main(int argc, char* argv[])
{
    // Defaults:
    int minDepth = 2;
    int maxDepth = 6;
    int iterations = 100;

    int worldSize = 4;

    circle c;

    double startPos[2] {0.3, 0.3};
    double endPos[2] {0.7, 0.7};

    c.r = 0.2;
    c.pos[0] = startPos[0];
    c.pos[1] = startPos[1];

    tftt::options.ghostsFlag = 0;
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
            tfetch("tftt.ghosts", tftt::options.ghostsFlag);
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
              << "\ttftt.ghosts = " << tftt::options.ghostsFlag << "\n"
              << "\ttftt.two2one = " << tftt::options.two2oneFlag << std::endl;

    // Init tree to min depth
    tftt::init(1.0, 1.0);

    for (int d = 0; d < minDepth; d++) {
        for (auto& cl : tftt::leaves) {
            tftt::refine(cl);
        }
    }

    tftt::drawMesh("circledata/mesh.min.dat");
    tftt::drawCurve("circledata/hilb.min.dat");
    tftt::drawBoundaries("circledata/bound.min.dat");

    // Refine to circle.
    for (int d = minDepth; d < maxDepth; d++) {
        tftt::adaptBegin();

        for (auto& cl : tftt::leaves) {
            if (c.intersects(cl))
                tftt::adaptAdd(cl);
        }

        tftt::adaptCommit();
    }

    tftt::drawMesh("circledata/mesh.init.dat");
    tftt::drawCurve("circledata/hilb.init.dat");
    tftt::drawBoundaries("circledata/bound.init.dat");


    double pos[2] = {0.37, 0.1};
    tftt::cell_t cl = tftt::atPos(pos);

    std::cout << cl << std::endl;

    std::ofstream sing("circledata/single.dat");
    tftt::drawCell(sing, cl);
    sing.close();

    tftt::calcFaceCoefs(cl);

    std::ofstream ngb("circledata/ngb.dat");
    tftt::drawPoissonNeighbourhood(ngb, cl);
    ngb.close();

    return 0;


    tftt::distribute(worldSize);
    tftt::saveParTree("circle.r{0}.tr", worldSize);

    // for (auto& cl : tftt::leaves) {
    //     if (cl.rank() == -1) {
    //         // Why?
    //         std::cout << cl << " = -1\n";
    //     }
    // }

    for (int n = 0; n < worldSize; n++) {
        std::cout << "Node = " << n << "\n";
        tftt::reset();
        tftt::loadParTree(formatString("circle.r{0}.tr", n));

        tftt::drawPartialMesh(formatString("parcircle/mesh.r{0}.dat", n));
        tftt::drawPartialCurve(formatString("parcircle/hilb.r{0}.dat", n));
        tftt::drawMesh(formatString("parcircle/mesh.r{0}.full.dat", n));
        tftt::drawCurve(formatString("parcircle/hilb.r{0}.full.dat", n));

        for (int b = 0; b < worldSize; b++) {
            tftt::drawGhosts(formatString("parcircle/ghosts.r{0}.b{1}.init.dat", n, b), b);
            tftt::drawBorder(formatString("parcircle/border.r{0}.b{1}.init.dat", n, b), b);
        }


        std::cout << "\tCells: " << tftt::gtree.ccells << "\n";
        std::cout << "\tGhosts: " << tftt::gtree.ghosts.size() << "\n";
    }

    int coarsened;
    bool cantCoarsen;
    tftt::cell_t tmp;
    for (ITER = 1; ITER <= iterations; ITER++) {
        c.pos[0] = startPos[0] + (endPos[0] - startPos[0])/iterations * ITER;
        c.pos[1] = startPos[1] + (endPos[1] - startPos[1])/iterations * ITER;

        std::cout << "Iteration " << ITER << "\n\tCircle at: " << c.pos[0] << "," << c.pos[1] << "\n";

        tftt::adaptBegin();
        for (auto& cl : tftt::leaves) {
            if (cl.level() < maxDepth && c.intersects(cl))
                tftt::adaptAdd(cl);
        }
        tftt::adaptCommit();

        std::cout << "\tAdapted " << tftt::adaptList.size() << "\n";

        do {
            coarsened = 0;
            tftt::adaptBegin();
            for (auto& cl : tftt::leaforthos) {
                if (c.intersects(cl))
                    continue;

                if (cl.level() > minDepth) {
                    cantCoarsen = false;
                    for (int ch = 0; ch < 1<<DIM; ch++) {
                        if (cl.child(ch).hasChildren()) {
                            cantCoarsen = true;
                            break;
                        }
                        cantCoarsen = tftt::findAround(cl, tftt::options.two2oneFlag-1,
                        [](tftt::cell_t& cl) {
                            return cl.hasGrandChildren();
                        });
                        if (cantCoarsen) break;
                    }

                    if (!cantCoarsen) {
                        tftt::adaptAddCoarsen(cl);
                    }
                }
            }
            tftt::adaptCommitCoarsen();

            std::cout << "\tCoarsened: " << tftt::adaptList.size() << "\n";
        }
        while (tftt::adaptList.size());

        std::cout << "\tCell Count: " << tftt::gtree.ccells << "\n";

        // int ccells = 0;
        // for (auto& cl : tftt::leaves) {
        //     ccells++;
        // }
        // std::cout << "\tCell Count (from leaves): " << ccells << "\n";

        // ccells = 0;
        // for (auto& cl : tftt::curve) {
        //     ccells++;
        // }
        // std::cout << "\tCell Count (from curve): " << ccells << "\n";

        tftt::drawMesh(formatString("circledata/mesh.{0}.dat", ITER));
        tftt::drawCurve(formatString("circledata/hilb.{0}.dat", ITER));
        tftt::drawBoundaries(formatString("circledata/bound.{0}.dat", ITER));
    }

    return 0;
}
