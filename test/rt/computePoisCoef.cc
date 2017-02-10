
#include "tftt/tftt.h"

#include "phys.h" // For rho1, rho2


double faceVoF(tftt::cell_t const& c1, int dir) {
    return 1.0; // Todo: VoF calculation.
}

double faceRho(tftt::cell_t const& c1, int dir) {
    double vof = faceVoF(c1, dir);
    return rho1*vof + (1.0-vof)*rho2;
}


void computePoisCoef()
{
    double dx, dy, rho, sumCoef;
    int i;

    // x,y coefficients
    // TODO: Re-use results. Could be done by storing time computed?
    for (auto& cl : tftt::leaves) {
        dx = cl.size(0);
        dy = cl.size(1);
        sumCoef = 0.0;

        // X coefs, both directions
        for (i = 0; i <=1; i++) {
            rho = faceRho(cl, i);
            // cl.facedata(i).poisCoef = dy*dy / (dx*dy * rho);
            sumCoef += // result of below
            cl.facedata(i).poisCoef = dy / (dx * rho);
        }

        // Y coefs, both directions
        for (i = 2; i <= 3; i++) {
            rho = faceRho(cl, i);
            // cl.facedata(i).poisCoef = dx*dx / (dx*dy * rho);
            sumCoef += // result of below
            cl.facedata(i).poisCoef = dx / (dy * rho);
        }

        // Central Coefficient
        if (sumCoef == 0.0)
            cl->cenCoef = 1.0;
        else
            cl->cenCoef = 1.0 / sumCoef;
    }
}

