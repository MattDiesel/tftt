
#include <iostream>
#include <fstream>
#include <algorithm>

#include "util/formatstring.h"

#include "tftt/tftt.h"


using namespace util; // for formatString


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

double dist(tftt::cell_t a, tftt::cell_t b)
{
    double x = a.centre(0) - b.centre(0);
    double y = a.centre(1) - b.centre(1);
    return x*x + y*y;
}


void drawCurves(std::string fnameFmt)
{
    int node = 0;
    std::ofstream ofs(formatString(fnameFmt, node));

    for (auto cl : tftt::curve) {
        if (cl.rank() != node) {
            node++;
            ofs.close();
            ofs.open(formatString(fnameFmt, node));
        }

        ofs << cl.centre(0) << " " << cl.centre(1) << "\n";
    }

    ofs.close();
}


int main(int argc, char* argv[])
{
    // Defaults:
    int minDepth = 3, maxDepth = 5;
    constexpr int worldSize = 3;
    int pass = 20;

    circle c;
    c.r = 0.2;
    c.pos[0] = 0.3;
    c.pos[1] = 0.3;

    tftt::options.two2oneFlag = 3;

    std::cout << "Using Parameters: \n"
              << "\tpass = " << pass << "\n"
              << "\tminDepth = " << minDepth << "\n"
              << "\tmaxDepth = " << maxDepth << "\n"
              << "\tworldSize = " << worldSize << "\n"
              << "\tc.pos[0] = " << c.pos[0] << "\n"
              << "\tc.pos[1] = " << c.pos[1] << "\n"
              << "\tc.r = " << c.r << "\n"
              << "\ttftt.two2one = " << tftt::options.two2oneFlag << std::endl;

    // Init tree to min depth
    tftt::init(1.0, 1.0);
    for (int d = 0; d < minDepth; d++) {
        for (auto& cl : tftt::leaves) {
            tftt::refine(cl);
        }
    }

    // Refine to circle.
    for (int d = minDepth; d < maxDepth; d++) {
        tftt::adaptSwBegin();
        for (auto& cl : tftt::leaves) {
            if (c.intersects(cl))
                tftt::adaptSwSetRefine(cl);
        }
        tftt::adaptSwCommit();
    }

    for (auto& cl : tftt::leaves) {
        tftt::calcFaceCoefs(cl);
    }

    tftt::distribute(worldSize);

    tftt::cell_t firstRank1;

    for (auto cl : tftt::curve) {
        if (cl.rank() == 1) {
            firstRank1 = cl;
            break;
        }
    }
    std::cout << "Changeover cell is at " << firstRank1 << "\n";

    tftt::cell_t lastSend = firstRank1.prev();
    for (int i = 0; i < pass; i++) {
        lastSend = lastSend.prev();
    }
    std::cout << "New rank 1 start at " << lastSend << "\n";


    // Fig 1: Initial curve
    tftt::drawPrettyMesh("report/movecells/mesh.init.dat");
    drawCurves("report/movecells/hilb.r{0}.init.dat");

    // Background heatmap hack
    {
        int steps = 2 << maxDepth;
        std::ofstream heatmap("report/ghosts/background.dat");
        tftt::cell_t cl;
        double scale = 1.0 / steps;
        double p[2];
        for (int y = 0; y < steps; y++) {
            for (int x = 0; x < steps; x++) {
                p[0] = x*scale+scale*0.5;
                p[1] = y*scale+scale*0.5;
                cl = tftt::atPos(p);
                heatmap << p[0] << " " << p[1] << " " << (int)cl.rank() << "\n";
            }
            heatmap << "\n";
        }
        heatmap.close();
    }


    // Fig 1.5: Ghost and border cells
    {
        tftt::drawPrettyMesh("report/ghosts/mesh.init.dat");

        std::set<tftt::cell_t, tftt::cell_t::parless> Gp[worldSize][worldSize];
        std::set<tftt::cell_t, tftt::cell_t::parless> B[worldSize];

        for (auto cl : tftt::curve) {
            for (auto P : cl.poissonNeighbourhood()) {
                if (P.isBoundary()) continue;

                if (P.rank() != cl.rank()) {
                    Gp[cl.rank()][P.rank()].insert(P);
                    B[cl.rank()].insert(cl);
                }
            }
        }

        {
            for (int r = 0; r < worldSize; r++) {
                for (int b = 0; b < worldSize; b++) {
                    if (r == b) continue;

                    std::ofstream gh(formatString("report/ghosts/gh.r{0}.b{1}.dat", r, b));
                    std::ofstream bd(formatString("report/ghosts/bd.r{1}.b{0}.dat", r, b));
                    for (auto g : Gp[r][b]) {
                        tftt::drawCell(gh, g);
                        tftt::drawCell(bd, g);
                    }
                    gh.close();
                    bd.close();
                }
            }
        }

        tftt::cell_t current, next;
        double curDist;
        for (int r = 0; r < worldSize; r++) {
            if (B[r].empty()) continue;

            std::ofstream bd(formatString("report/ghosts/bd.r{0}.dat", r));

            // Start at the leftmost
            current = tftt::cell_t();
            for (auto g : B[r]) {
                if (!current.isValid()) {
                    current = g;
                }
                else if (r == 2) {
                    if (g.centre(1) > current.centre(1))
                        current = g;
                }
                else if (g.centre(0) == current.centre(0)) {
                    if (g.centre(1) > current.centre(1))
                        current = g;
                }
                else if (g.centre(0) < current.centre(0)) {
                    current = g;
                }
            }

            B[r].erase(current);
            bd << current.centre(0) << " " << current.centre(1) << "\n";

            while (!B[r].empty()) {
                curDist = 20.0;
                for (auto g : B[r]) {
                    if (dist(current, g) < curDist) {
                        next = g;
                        curDist = dist(current, g);
                    }
                }
                current = next;

                B[r].erase(current);
                bd << current.centre(0) << " " << current.centre(1) << "\n";
            }

            bd.close();
        }
    }

    // Fig 2: Exchange
    tftt::plot::partialHilbert("report/movecells/hilb.overlapl.dat",
                           lastSend, firstRank1.prev());
    tftt::plot::partialHilbert("report/movecells/hilb.overlapr.dat",
                           lastSend.next(), firstRank1);

    {
        std::ofstream ghostNotify("report/movecells/ghostnotify.dat");
        tftt::cell_t cl = lastSend;
        do {
            for (auto nb : cl.neighbours()) {
                if (nb.rank() > 1) {
                    tftt::drawCell(ghostNotify, cl);
                    break;
                }
            }
            cl = cl.next();
        }
        while (cl != firstRank1);
    }

    tftt::cell_t cl = firstRank1.prev();
    for (int i = 0; i < pass; i++) {
        cl.rank() = 1;
        cl = cl.prev();
    }

    {
        std::ofstream Tr("report/movecells/tr.dat");

        cl = firstRank1.prev();
        for (int i = 0; i < pass; i++) {
            tftt::drawCell(Tr, cl);

            tftt::TreeCell& tc = cl.group->cells[cl.index];

            for (int p = 0; p < tc.poisNgbC; p++) {
                if (tc.poisNgb[p].rank() != 1) {
                    tftt::drawCell(Tr, tc.poisNgb[p]);
                }
            }

            cl = cl.prev();
        }
    }

    // Fig N: Final
    drawCurves("report/movecells/hilb.r{0}.final.dat");



    return 0;
}
