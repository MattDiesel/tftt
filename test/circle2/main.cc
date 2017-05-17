
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

    tftt::drawMesh("circle2data/mesh.min.dat");

    // Refine to circle.
    for (int d = minDepth; d < maxDepth; d++) {
        tftt::adaptSwBegin();

        for (auto& cl : tftt::leaves) {
            if (c.intersects(cl))
                tftt::adaptSwSetRefine(cl);
        }

        tftt::adaptSwCommit();
    }

    // tftt::adaptSwBegin();
    // double p[2] = {0.5, 0.5};
    // double p2[2] = {0.6, 0.6};
    // tftt::adaptSwSetRefine(tftt::atPos(p));
    // tftt::adaptSwCommit();

    // tftt::adaptSwBegin();
    // tftt::adaptSwSetRefine(tftt::atPos(p));
    // tftt::adaptSwCommit();

    // tftt::adaptSwBegin();
    // tftt::adaptSwSetRefine(tftt::atPos(p2));
    // tftt::adaptSwCommit();

    // tftt::adaptSwBegin();
    // tftt::adaptSwSetRefine(tftt::atPos(p2));
    // tftt::adaptSwCommit();

    tftt::drawMesh("circle2data/mesh.init.dat");
    // tftt::drawCurve("circle2data/hilb.init.dat");
    // tftt::drawBoundaries("circle2data/bound.init.dat");

    // return 0;

    tftt::cell_t tmp;
    for (ITER = 1; ITER <= iterations; ITER++) {
        c.pos[0] = startPos[0] + (endPos[0] - startPos[0])/iterations * ITER;
        c.pos[1] = startPos[1] + (endPos[1] - startPos[1])/iterations * ITER;

        std::cout << "Iteration " << ITER << "\n\tCircle at: " << c.pos[0] << "," << c.pos[1] << "\n";

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
            if (cl.level() > minDepth) {
                tftt::adaptSwSetCoarsen(cl);
            }
        }
        tftt::adaptSwCommit();

        // do {
        //     coarsened = 0;
        //     tftt::adaptBegin();
        //     for (auto& cl : tftt::leaforthos) {
        //         if (c.intersects(cl))
        //             continue;

        //         if (cl.level() > minDepth) {
        //             cantCoarsen = false;
        //             for (int ch = 0; ch < 1<<DIM; ch++) {
        //                 if (cl.child(ch).hasChildren()) {
        //                     cantCoarsen = true;
        //                     break;
        //                 }
        //                 cantCoarsen = tftt::findAround(cl, tftt::options.two2oneFlag-1,
        //                 [](tftt::cell_t& cl) {
        //                     return cl.hasGrandChildren();
        //                 });
        //                 if (cantCoarsen) break;
        //             }

        //             if (!cantCoarsen) {
        //                 tftt::adaptAddCoarsen(cl);
        //             }
        //         }
        //     }
        //     tftt::adaptCommitCoarsen();

        //     std::cout << "\tCoarsened: " << tftt::adaptList.size() << "\n";
        // }
        // while (tftt::adaptList.size());

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

        tftt::drawMesh(formatString("circle2data/mesh.{0}.dat", ITER));
        tftt::drawCurve(formatString("circle2data/hilb.{0}.dat", ITER));
        tftt::drawBoundaries(formatString("circle2data/bound.{0}.dat", ITER));
    }

    return 0;
}
