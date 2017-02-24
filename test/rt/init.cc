
#include <cmath>
#include <iostream>

#include "tftt/tftt.h"
#include "util/pars.h"

#include "boundarycond.h"
#include "readData.h"


using namespace util; // tfetch


void initCondition() {
    const double pi = 3.1415926537;
    double x, y, high, ampl, wavenr, alpha; 
    int i;


    tfetch("amplitude", ampl);
    tfetch("wavenr", wavenr);
    alpha = 2.0*wavenr*pi;

    std::cout << "Initial condition for Rayleigh-Taylor Instability" << std::endl;
    std::cout << "amplitude = " << ampl << " alpha = " << alpha << std::endl;

    // for(auto& cl : tftt::leaves) {
    //     cl->u = 0.0;
    //     cl->v = 0.0;
    // }

    // Split is a sine curve at 2.0
    // Refine to geometry, and calculate VoF. 


    // Refine to mininum mesh size. Default = 2
    int mindepth = 2, maxdepth = 8;
    tfetch("mindepth", mindepth);
    tfetch("mindepth", maxdepth);

    for (i = 0; i < mindepth; i++) {
        for (auto& cl : tftt::leaves) {
            tftt::refine(cl);
        }
    }

    // Refine cells straddling the curve.
    double yl, yr;
    for (i = mindepth; i < maxdepth; i++) {
        tftt::adaptBegin();
        for (auto& cl : tftt::leaves) {
            yl = ampl*cos(alpha*cl.vertex(0, 0))+2.0;
            yr = ampl*cos(alpha*cl.vertex(1, 0))+2.0;

            if ((cl.vertex(0, 1) < yl && cl.vertex(2, 1) > yl)
                    || (cl.vertex(1, 1) < yr && cl.vertex(3, 1) > yr)) {
                tftt::adaptAdd(cl);
            }
        }
        tftt::adaptCommit();
    }

    tftt::drawMesh("init.mesh");

    // Set VoF for all cells

    // NB: Not 100% accurate, just averages cell sides. This only differs in the
    //     case that the curve passes through the top/bottom, but these cells 
    //     should be small.
    double fl, fr;
    for (i = mindepth; i < maxdepth; i++) {
        for (auto& cl : tftt::leaves) {
            cl->cc = 0.0;

            for (int d = 0; d < 2; d++) {
                y = ampl*cos(alpha*cl.vertex(d, 0))+2.0;
                y = cl.vertex(2, 1) - y;
                y /= cl.size(1);

                cl->cc += std::min(std::max(y, 0.0), 1.0);
            }
            cl->cc *= 0.5;

            cl->v[0] = 0.0;
            cl->v[1] = 0.0;
        }
    }

    tftt::drawMatrix("init.matrix.pgm", 256, 1024, [](tftt::data_t& dt, int max) {
        return dt.cc*max;
    });

    bcV();
    bcVof();
}

