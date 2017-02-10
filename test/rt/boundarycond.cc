
#include "tftt/tftt.h"

#include "readData.h"


/* boundary conditions for velocity */
// periodic condition in x-direction and slip condition on top, bottom wall
void bcV() {

    // No Slip on top/bottom
    for (int b = 0; b < 2*DIM; b++) {
        for (auto& cl : tftt::boundaryCells(b)) {
            cl->v[0] = cl.neighbour(b ^ 1)->v[0];
            cl->v[1] = 0.0;
        }
    }

    // Think there should be wrap-around on left/right. Implementation tricky.
}

void bcPressure() {
    for (int b = 0; b < 2*DIM; b++) {
        for (auto& cl : tftt::boundaryCells(b)) {
            cl->P = cl.neighbour(b ^ 1)->P;
        }
    }

    // Again I think there's supposed to be some sort of wrap around.
}

void bcVof() {
    
    // Constrain
    for (auto& cl : tftt::leaves) {
        cl->cc = std::min(std::max(cl->cc, 0.0), 1.0);
    }

    // Copy boundaries
    for (int b = 0; b < 2*DIM; b++) {
        for (auto&cl : tftt::boundaryCells(b)) {
            cl->cc = cl.neighbour(b ^ 1)->cc;
        }
    }
}
