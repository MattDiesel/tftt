
#include <iostream>
#include <fstream>

#include "util/formatstring.h"
#include "util/pars.h"

#include "tftt/tftt.h"
#include "tftt/tree.h"

using namespace util; // for formatString


struct circle {

    double pos[2];
    double r;

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
    // Defaults:
    int minDepth = 2;
    int maxDepth = 6;
    int iterations = 100;
    int plotEvery = 1;

    double dirichlet = 0.0;
    double neuman = 0.0;
    tftt::options.isNeuman = true;

    double omega = 1.3;

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

            tfetch("plotEvery", plotEvery);
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

            if (tfetch("neuman", neuman))
                tftt::options.isNeuman = true;
            if (tfetch("dirichlet", dirichlet))
                tftt::options.isNeuman = false;
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

    // tftt::drawMesh("pois/mesh.min.dat");
    // tftt::drawCurve("pois/hilb.min.dat");
    // tftt::drawBoundaries("pois/bound.min.dat");

    // Refine to circle.
    for (int d = minDepth; d < maxDepth; d++) {
        tftt::adaptBegin();

        for (auto& cl : tftt::leaves) {
            if (c.intersects(cl))
                tftt::adaptAdd(cl);
        }

        tftt::adaptCommit();
    }

    // Initial Conditions for cells
    double fl, fr;
    for (auto& cl : tftt::leaves) {
        cl->cc = c.cc(cl);
        cl->P = dirichlet;
        tftt::calcFaceCoefs(cl);
    }

    if (!tftt::options.isNeuman) {
        for (int b = 0; b < 4; b++) {
            for (auto& cl : tftt::boundaryCells(b)) {
                cl->P = dirichlet;
            }
        }
    }

    tftt::drawMesh("pois/mesh.init.dat");
    tftt::drawCurve("pois/hilb.init.dat");
    tftt::drawBoundaries("pois/bound.init.dat");


    // tftt::drawMatrix("pois/cc.init.pgm", 512, 512, [](tftt::data_t& dt, int max) {
    //     return dt.cc*max;
    // });


    // double pos[2] = {0.148438,0.425};
    // tftt::cell_t cl = tftt::atPos(pos);

    // std::cout << cl << " is " << cl.id().id << std::endl;

    // std::ofstream sing("pois/single.dat");
    // tftt::drawCell(sing, cl);
    // sing.close();

    // tftt::calcFaceCoefs(cl);

    // std::ofstream ngb("pois/ngb.dat");
    // tftt::drawPoissonNeighbourhood(ngb, cl);
    // ngb.close();

    // tftt::drawMatrix("pois/cc.init.pgm", 512, 512, [](tftt::data_t& dt, int max) {
    //     return dt.cc*max;
    // });
    // tftt::drawMatrix("pois/P.init.pgm", 512, 512, [](tftt::data_t& dt, int max) {
    //     return (dt.P+1.0)*max * 0.5;
    // });


    // {
    //     char c;
    //     std::cin >> c;
    // }

    std::ofstream ofGraph("pois/res.graph.dat");

    tftt::cell_t mx;
    for (ITER = 0; ITER < iterations; ITER++) {
        std::cerr << "Iteration " << ITER << "\n";
        mx = tftt::max(dtPc);
        std::cerr << "Max P = " << mx->P << " @ " << mx << "\n";

        mx = tftt::max(dtPcn);
        std::cerr << "Min P = " << mx->P << " @ " << mx << "\n";

        tftt::relax(omega, dtP, fn);

        if (!tftt::options.isNeuman) {
            for (int b = 0; b < 4; b++) {
                for (auto& cl : tftt::boundaryCells(b)) {
                    cl->P = dirichlet;
                }
            }
        }

        ofGraph << tftt::resid(dtP, fn) << "\n";

        std::cout << std::endl;

        if (ITER % plotEvery == 0) {
            tftt::plotMatrix(formatString("pois/P.{0}.dat", ITER), [](tftt::data_t& dt) {
                return dt.P;
            });

            tftt::plotMatrix(formatString("pois/res.{0}.dat", ITER), [](tftt::data_t& dt) {
                return dt.res;
            });
        }
    }

    tftt::plotMatrix("pois/P.final.dat", [](tftt::data_t& dt) {
        return dt.P;
    });
    tftt::plotMatrix("pois/res.final.dat", [](tftt::data_t& dt) {
        return dt.res;
    });

    return 0;
}
