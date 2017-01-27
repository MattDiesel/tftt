
#define TFTT_DEBUG


#include <iostream>
#include <set>

#include "tftt/tftt.h"

#include "tftt/treegroup.h"
#include "formatstring.h"


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


double faceFlux(tftt::cell_t cl, tftt::cell_t nb) {

    double nbDat = 0.0;
    if (nb.hasChildren()) {
        for (auto& ch : *nb.children()) {
            nbDat += ch.data();

            if (ch.data() == 1337) {
                throw;
            }
        }
        nbDat *= 0.25;
    }
    else {
        nbDat = nb.data();

        if (nbDat == 1337) {
            throw;
        }
    }

    return cl.data() - nbDat;
}




int main(int argc, char const *argv[])
{
	std::cout << "Init Tree" << std::endl;

	tftt::init(4.0, 2.0);

    std::cout << "Min Depth" << std::endl;

    int minDepth = 3;
    for (int i = 0; i < minDepth; i++) {
        for (auto& cl : tftt::leaves) {
            refine(cl);
        }
    }

    tftt::drawMesh("minDepth.dat");


    std::cout << "Refinement" << std::endl;

    for (auto& cl : tftt::leaves) {
        cl.data() = (mbrot(cl.centre(0) - 2.0, cl.centre(1) - 1.0));
    }


    int iter = 0;
    do {
        tftt::adaptBegin();

        int maxDepth = 10;
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
                ch.data() = (mbrot(ch.centre(0) - 2.0, ch.centre(1) - 1.0));
            }
        }

        iter++;

        std::cout << "Iteration "
            << iter << " refined: "
            << tftt::adaptList.size() << std::endl;

    } while (!tftt::adaptList.empty());

    tftt::drawMesh(formatString("refined{0}.dat", iter));
    tftt::drawMatrix(formatString("refined{0}.pgm", iter), 1024, 512);




	return 0;
}
