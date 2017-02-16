
#include <iostream>

#include "formatstring.h"
#include "tftt/tftt.h"
#include "tftt/tree.h"

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
            if (contains(cl.vertex(v, 0), cl.vertex(v, 1))) {
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



int main(int argc, char const *argv[])
{

    int minDepth = 2;
    int maxDepth = 6;
    int iterations = 100;

    circle c;

    double startPos[2] {0.3, 0.3};
    double endPos[2] {0.7, 0.7};

    c.r = 0.2;
    c.pos[0] = startPos[0];
    c.pos[1] = startPos[1];

    // Init tree to min depth
    tftt::init(1.0, 1.0);

    for (int d = 0; d < minDepth; d++) {
        for (auto& cl : tftt::leaves) {
            tftt::refine(cl);
        }
    }

    tftt::drawMesh("circledata/mesh.min.dat");
    tftt::drawCurve("circledata/hilb.min.dat");
    tftt::drawBoundaries("circledata/bound.min.dat");

    // Refine to circle. 
    for (int d = minDepth; d < maxDepth; d++) {
        tftt::adaptBegin();

        for (auto& cl : tftt::leaves) {
            if (c.intersects(cl))
                tftt::adaptAdd(cl);
        }

        tftt::adaptCommit();
    }

    tftt::drawMesh("circledata/mesh.init.dat");
    tftt::drawCurve("circledata/hilb.init.dat");
    tftt::drawBoundaries("circledata/bound.init.dat");

    int coarsened;
    bool cantCoarsen;
    tftt::cell_t tmp;
    for (ITER = 1; ITER <= iterations; ITER++) {
        c.pos[0] = startPos[0] + (endPos[0] - startPos[0])/iterations * ITER;
        c.pos[1] = startPos[1] + (endPos[1] - startPos[1])/iterations * ITER;

        std::cout << "Iteration " << ITER << "\n\tCircle at: " << c.pos[0] << "," << c.pos[1] << "\n";

        tftt::adaptBegin();
        for (auto& cl : tftt::leaves) {
            if (cl.level() < maxDepth && c.intersects(cl))
                tftt::adaptAdd(cl);
        }
        tftt::adaptCommit();

        std::cout << "\tAdapted " << tftt::adaptList.size() << "\n";

        coarsened = 0;
        tftt::adaptBegin();
        for (auto& cl : tftt::leaforthos) {
            if (c.intersects(cl))
                continue;

            cantCoarsen = false;
            for (int ch = 0; ch < 1<<DIM; ch++) {
                if (cl.child(ch).hasChildren()) {
                    cantCoarsen = true;
                    break;
                }
                for (int nb = 0; nb < 2*DIM; nb++) {
                    tmp = cl.child(ch).neighbour(nb);
                    if (tmp.hasChildren()) {
                        cantCoarsen = true;
                        break;
                    }
                }
                if (cantCoarsen) break;
            }

            if (!cantCoarsen) {
                tftt::adaptAddCoarsen(cl);
            }
        }
        tftt::adaptCommitCoarsen();

        std::cout << "\tCoarsened: " << tftt::adaptList.size() << "\n";

        std::cout << "\tCell Count: " << tftt::gtree.ccells << "\n";

        tftt::drawMesh(formatString("circledata/mesh.{0}.dat", ITER));
        tftt::drawCurve(formatString("circledata/hilb.{0}.dat", ITER));
        tftt::drawBoundaries(formatString("circledata/bound.{0}.dat", ITER));
    }

    return 0;
}
