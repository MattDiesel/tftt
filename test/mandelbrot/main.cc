
#define TFTT_DEBUG


#include <iostream>
#include <set>

#include "tftt/tftt.h"
#include "tftt/tree.h"
#include "tftt/treegroup.h"

#include "formatstring.h"
#include "pars.h"


int escape_iterations(double re, double im, int max_iterations) {
    double xre = 0.0, xim = 0.0, rtemp, xre2 = 0.0, xim2 = 0.0;
    int i = max_iterations;

    do {
        --i;
        rtemp = xre2 - xim2 + re;
        xim = 2*xre*xim + im;
        xre = rtemp;
        xre2 = rtemp*rtemp;
        xim2 = xim*xim;
    } while ( (xre2 + xim2 < 4.0f) && i);

    return max_iterations-i;
}

double mbrot(double re, double im) {
    return escape_iterations(re, im, 255) / 255.0;
}


double avrChildren(tftt::cell_t cl) {
    double ret = 0.0;
    for (auto& ch : *cl.children()) {
        if (ch.hasChildren()) {
            ret += avrChildren(ch);
        }
        else {
            ret += ch->P;
        }
    }

    return ret*0.25;
}

double faceFlux(tftt::cell_t cl, tftt::cell_t nb) {

    double nbDat = 0.0;
    if (nb.hasChildren()) {
        nbDat = avrChildren(nb);
    }
    else {
        nbDat = nb->P;

        if (nbDat == 1337) {
            throw;
        }
    }

    return cl->P - nbDat;
}




int main(int argc, char const *argv[])
{
    int cnodes = 4;
    int minDepth = 3;
    int maxDepth = 8;

    if (argc > 1) {
        std::cout << "Reading parameters from: " << argv[1] << "\n";
        getpars(argv[1]);
        tfetch("minDepth", minDepth);
        tfetch("maxDepth", maxDepth);
        tfetch("cnodes", cnodes);
    }

    std::cout << "Using configuration:\n"
            << "\tminDepth = " << minDepth << "\n"
            << "\tmaxDepth = " << maxDepth << "\n"
            << "\tcnodes = " << cnodes << std::endl;


	std::cout << "Init Tree" << std::endl;

	tftt::init(4.0, 2.0);

    std::cout << "Min Depth" << std::endl;

    for (int i = 0; i < minDepth; i++) {
        for (auto& cl : tftt::leaves) {
            refine(cl);
        }
    }

    // tftt::drawMesh("minDepth.dat");

    std::cout << "Refinement" << std::endl;

    for (auto& cl : tftt::leaves) {
        cl->P = (mbrot(cl.centre(0) - 2.0, cl.centre(1) - 1.0));
    }

    int iter = 0;
    do {
        tftt::adaptBegin();

        double eps;
        tftt::cell_t nb;
        for (auto& cl : tftt::leaves) {
            if (cl.level() >= maxDepth) 
                continue;

            for (int n = 0; n < 2*DIM; n++) {
                nb = cl.neighbour(n);

                if (nb.isBoundary()) {
                    continue;
                }

                eps = faceFlux(cl, nb);
                if (eps < -0.1 || eps > 0.1) {
                    tftt::adaptAdd(cl);
                }
            }
        }

        tftt::adaptCommit();

        // Reprocess new cells
        for (auto& cl : tftt::adaptList) {
            for (auto& ch : *cl.children()) {
                ch->P = (mbrot(ch.centre(0) - 2.0, ch.centre(1) - 1.0));
            }
        }

        iter++;

        std::cout << "Iteration "
            << iter << " refined: "
            << tftt::adaptList.size() << std::endl;

    } while (!tftt::adaptList.empty());

    std::cout << "Refined to depth " << maxDepth << "\n";
    std::cout << "Total Cells: " << tftt::gtree.ccells << "\n";

    tftt::drawMesh("mbrot/mesh.init.dat");
    tftt::drawCurve("mbrot/hilb.init.dat");
    tftt::drawBoundaries("mbrot/bound.init.dat");
    tftt::drawMatrix("mbrot/init.pgm", 1024, 512, [](tftt::data_t& dt, int max) {
        return dt.P;
    });

    tftt::distribute(cnodes);
    tftt::splitToDisk("mbrot/r{0}.tr");

    for (int n = 0; n < cnodes; n++) {
        std::cout << "Node = " << n << "\n";
        tftt::reset();
        tftt::loadTree(formatString("mbrot/r{0}.tr", n), n);

        tftt::drawPartialMesh(formatString("mbrot/mesh.r{0}.dat", n));
        tftt::drawPartialCurve(formatString("mbrot/hilb.r{0}.dat", n));
        tftt::drawGhosts(formatString("mbrot/ghosts.r{0}.dat", n));
        tftt::drawBoundaries(formatString("mbrot/bound.r{0}.dat", n));

        std::cout << "\tGhosts: " << tftt::gtree.ghosts.size() << "\n";
    }

	return 0;
}
