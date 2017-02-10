
#include "tftt/tftt.h"

#include "phys.h"
#include "readData.h"
#include "boundarycond.h"

// Generate the lambdas easily
#define DTFN(x) (tftt::fnData) [](tftt::data_t& dt) -> double { return dt. x; }


double interpNegDiag(tftt::cell_t cl, tftt::fnData dt) {
    double ret;

    if (cl.neighbour(0).isBoundary() && cl.neighbour(2).isBoundary()) {
        return dt(cl.data());
    }

    // No simple way to interpolate diagonal cell atm.
    bool lless = cl.neighbour(0).level() < cl.level();
    bool bless = cl.neighbour(2).level() < cl.level();
    if (lless && bless) { // Both left and below less refined
        tftt::cell_t diag = cl.neighbour(0).neighbour(2);
        if (diag.hasChildren())
            ret = dt(diag.child(3).data());
        else {
            // This is really awful. Basically nothing to go on. 
            // Instead interpolate on the larger grid
            ret = cl.parent().avrChildren(dt);
            ret += dt(cl.neighbour(0).data());
            ret += dt(cl.neighbour(2).data());
            ret += dt(cl.neighbour(2).neighbour(0).data());
            ret *= 0.25;
        }
    }
    else if (lless) { // Use below neigbour, then left
        ret = cl.neighbour(2).ngbVal(0, dt);
    }
    else { // Use left, then below
        ret = cl.neighbour(0).ngbVal(2, dt);
    }

    return ret;
}




void computeRho() {
    double vof;

    for (auto& cl : tftt::leaves) {
        vof = 0.0;

        if (cl.neighbour(0).isBoundary() && cl.neighbour(2).isBoundary()) {
            // This case isn't possible with the boundary layer.
            vof = cl->cc;
        }
        else{
            vof += cl->cc;
            vof += cl.ngbVal(0, DTFN(cc));
            vof += cl.ngbVal(2, DTFN(cc));
            vof += interpNegDiag(cl, DTFN(cc));
            vof *= 0.25;
        }

        cl->rho = vof*rho1 + (1.0-vof)*rho2;
    }
}


void computeAdv() {
    tftt::fnData dtv[2] = {
        [](tftt::data_t& dt) -> double { return dt.v[0]; },
        [](tftt::data_t& dt) -> double { return dt.v[1]; }
    };

    for (int d = 0; d < DIM; d++) {
        for (auto& cl : tftt::leaves) {
            cl->F1[d] = 0.5 * (
                    cl->v[0] * ( cl.ngbVal(0, dtv[d]) - cl.ngbVal(1, dtv[d]) ) / cl.size(0)
                    + cl->v[1] * ( cl.ngbVal(2, dtv[d]) - cl.ngbVal(3, dtv[d]) ) / cl.size(1)
                );
        }
    }
}


void computeVisc() {
    double vof, mu;

    for (int d = 0; d < DIM; d++) {
        for (auto& cl : tftt::leaves) {
            vof = 0.5*( cl->cc + cl.ngbVal(2*d, DTFN(cc)) );
            mu = mu1*vof+(1.0-vof)*mu2;

            cl->D1[d] = mu*( cl.ngbVal(2*d | 1, DTFN(v[0])) - cl->v[0] ) / cl.size(d);
            cl->D2[d] = mu*( cl.ngbVal(2*d | 1, DTFN(v[1])) - cl->v[1] ) / cl.size(d);
        }
    }
}


void computeIntermVelo() {
    double temp;

    tftt::fnData dtD1[2] = {
        [](tftt::data_t& dt) -> double { return dt.D1[0]; },
        [](tftt::data_t& dt) -> double { return dt.D1[1]; }
    };
    tftt::fnData dtD2[2] = {
        [](tftt::data_t& dt) -> double { return dt.D2[0]; },
        [](tftt::data_t& dt) -> double { return dt.D2[1]; }
    };

    for (auto& cl : tftt::leaves) {
        for (int d = 0; d < DIM; d++) {
            temp = 0.0;

            // Todo: Interpolation if neighbours are at different levels
            temp += (cl->D1[d] - cl.ngbVal(0, dtD1[d])) / cl.size(0);
            temp += (cl->D2[d] - cl.ngbVal(2, dtD2[d])) / cl.size(1);

            temp /= cl->rho;
            temp += cl->F1[d];
            temp += g[d];
            temp *= tau;

            cl->v[d] = temp;
        }
    }

    // Reinforce boundary condition
    bcV();
}

void computeFlux() {
    double vl;

    tftt::fnData dtv[2] = {
        [](tftt::data_t& dt) -> double { return dt.v[0]; },
        [](tftt::data_t& dt) -> double { return dt.v[1]; }
    };

    for (int d = 0; d < DIM; d++) {
        for (auto& cl : tftt::leaves) {
            vl = tftt::interpALEVertex(cl, 1 << d, dtv[d]);
            cl->V[d] = 0.5*(cl->v[d]+vl) * cl.size(d ^ 1); // Todo: Generalise
        }
    }
}


//! compute divergence
//! Provides RHS of the Poisson Equation for pressure
void computeDive() {
    tftt::cell_t ngb;
    double u, v;

    tftt::fnData dtV[2] = {
        [](tftt::data_t& dt) -> double { return dt.V[0]; },
        [](tftt::data_t& dt) -> double { return dt.V[1]; }
    };

    for (auto& cl : tftt::leaves) {
        cl->dive = 0.0;

        for (int d = 0; d < DIM; d++) {
            ngb = cl.neighbour(d*2 | 1);
            if (ngb.isBoundary()) {
                // Boundary condition?
            }
            else {
                cl->dive += tftt::interpALEVertex(cl, 1 << d, dtV[d]);
            }
            cl->dive -= cl->V[d];
        }
    }
}

double computeCompatability() {
    double eps = 0.0;
    for (auto& cl : tftt::leaves) {
        eps += cl->dive;
    }
    return eps;
}


void correctVelo() {
    tftt::fnData dtP = [](tftt::data_t& dt) { return dt.P; };

    for (auto& cl : tftt::leaves) {
        cl->v[0] -= 0.5*(cl.ngbVal(2, dtP) - interpNegDiag(cl, dtP) - cl.ngbVal(0, dtP) + cl->P)
                    / (cl.size(0) * cl->rho);
        cl->v[1] -= 0.5*(cl.ngbVal(0, dtP) - interpNegDiag(cl, dtP) - cl.ngbVal(2, dtP) + cl->P)
                    / (cl.size(1) * cl->rho);
    }

    bcV();
}


void correctFlux() {
    tftt::fnData dtP = [](tftt::data_t& dt) { return dt.P; };

    for (auto& cl : tftt::leaves) {
        cl->V[0] += cl.facedata(0).poisCoef * (cl.ngbVal(0, dtP) - cl->P);
    }
    for (auto& cl : tftt::leaves) {
        cl->V[1] += cl.facedata(2).poisCoef * (cl.ngbVal(2, dtP) - cl->P);
    }
}
