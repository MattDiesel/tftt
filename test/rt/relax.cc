
#include "tftt/tftt.h"

#include "phys.h" // omega


void relax() {
    double temp, t1, t2;
    tftt::cell_t c1, c2;

    // Todo: Generalise.
    static const int c1map[4] = {1,0,2,0};
    static const int c2map[4] = {3,2,3,1};

    for (auto& cl : tftt::leaves) {
        temp = 0.0;

        for (int nb = 0; nb < 2*DIM; nb++) {
            // Todo: Boundary conditions currently enforced here. 
            //       Possibly more efficient to store boundary cells.
            if (cl.neighbour(nb).isBoundary()) {
                // Neumann condition, d/dx = 0
                temp += cl.facedata(nb).poisCoef * cl->P;
            }
            else if (cl.neighbour(nb).hasChildren()) {
                // Neighbour more refined
                // Use Neighbour cell's coefs on that border
                c1 = cl.neighbour(nb).child(c1map[nb]);
                c2 = cl.neighbour(nb).child(c2map[nb]);

                t1 = c1.facedata(nb ^ 1).poisCoef * c1->P;
                t2 = c2.facedata(nb ^ 1).poisCoef * c2->P;

                temp += (t1+t2) * 0.5;
            }
            else if (cl.neighbour(nb).level() < cl.level()) {
                // Neighbour less refined.
                // Interpolate for P, use coef for this cell
                temp += cl.facedata(nb).poisCoef 
                        * tftt::interpChild(cl.neighbour(nb), 
                                cl.index ^ (1 << (nb >> 1)),
                                nb ^ 1, [](tftt::data_t& dt) {
                                    return dt.P;
                                });
            }
            else {
                temp += cl.facedata(nb).poisCoef * cl.neighbour(nb)->P;
            }
        }

        temp -= cl->rhs;
        temp *= cl->cenCoef;

        cl->P = omega*temp + (1.0-omega)*cl->P;
    }
}
