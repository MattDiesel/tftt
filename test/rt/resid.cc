
#include "tftt/tftt.h"

#include "boundarycond.h"


void resid() {
    tftt::fnData dtP = [](tftt::data_t& dt) { return dt.P; };

    for (auto& cl : tftt::leaves) {
        cl->res = cl->rhs + cl->P / cl->cenCoef;

        for (int nb = 0; nb < 2*DIM; nb++) {
            if (cl.neighbour(nb).isBoundary()) {
                // Neuman d/dx = 0
                cl->res -= cl.facedata(nb).poisCoef * cl->P;
            }
            else {
                cl->res -= cl.facedata(nb).poisCoef * cl.ngbVal(nb, dtP);
            }
        }
    }
}


void rstrct() {
    tftt::fnData dtRes = [](tftt::data_t& dt) { return dt.res; };

    for (auto& cl : tftt::leaves) {
        cl->rhs = 0.25*( cl->res + cl.ngbVal(1, dtRes)
                + cl.ngbVal(3, dtRes) + cl.ngbVal(3, dtRes));
    }

    bcPressure();
}
