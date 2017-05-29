
#include <iostream>
#include <fstream>
#include <algorithm>

#include "util/formatstring.h"

#include "tftt/tftt.h"
#include "tftt/io/tikz.h"


using namespace util; // for formatString



// y = 0.3 + 6.189086*x - 6.682936*x^2 + 2.481516*x^3 - 0.2852701*x^4

double evalLine(double x)
{
    double coefs[5] = {0.3, 6.189086, -6.682936, 2.481516, -0.2852701};
    double xp = 1.0;

    double y = 0.0;
    for (int i = 0; i < 5; i++) {
        y += coefs[i]*xp;
        xp *= x;
    }

    return y;
}

void drawLayer(std::ostream& os, tftt::cell_t cl, int layer)
{
    if (cl.level() == layer) {
        tftt::drawCell(os, cl);
    }
    else if (cl.hasChildren()) {
        for (auto ch : *cl.children()) {
            drawLayer(os, ch, layer);
        }
    }
}

void drawLayer(std::string fname, int layer)
{
    std::ofstream ofs(fname);
    tftt::cell_t root = tftt::cell_t(-1);

    drawLayer(ofs, root, layer);
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


void tikzCurves(std::string fnameFmt)
{
    int node = 0;
    std::ofstream ofs(formatString(fnameFmt, node));

    node = -1;
    for (auto cl : tftt::curve) {
        if (node == -1) {
            node = 0;
            continue;
        }

        if (cl.rank() != node) {
            node++;
            ofs.close();
            ofs.open(formatString(fnameFmt, node));
        }
        else {
            ofs << "\\draw (" << cl.prev().centre(0) << "," << cl.prev().centre(1)
                << ") -- (" << cl.centre(0) << "," << cl.centre(1) << ");\n";
        }
    }

    ofs.close();
}


int main(int argc, char* argv[])
{
    // Defaults:
    int minDepth = 0, maxDepth = 5;
    int worldSize = 3;

    double p[2];

    std::string outDir = "../docs/report/method/gen";


    {
        std::ofstream ofline("report/ftt/line.dat");

        for (double x = 0.01; x < 3.79; x += 0.01) {
            p[0] = x;
            p[1] = evalLine(x);
            ofline << p[0] << " " << p[1] << "\n";
        }
    }

    for (int two2one = 0; two2one < 4; two2one++) {
        tftt::options.two2oneFlag = two2one;

        if (two2one > 0) {
            tftt::reset();
        }

        // Init tree to min depth
        tftt::init(4.0, 4.0);
        for (int d = 0; d < minDepth; d++) {
            for (auto& cl : tftt::leaves) {
                tftt::refine(cl);
            }
        }

        tftt::cell_t cl;
        for (int d = minDepth; d < maxDepth; d++) {
            tftt::adaptSwBegin();
            for (double x = 0.01; x < 3.79; x += 0.01) {
                p[0] = x;
                p[1] = evalLine(x);

                cl = tftt::atPos(p);
                if (cl.isValid()) {
                    tftt::adaptSwSetRefine(cl);
                }
            }
            tftt::adaptSwCommit();
        }

        tftt::tikz::mesh(formatString("{0}/mesh-t{1}.tex", outDir, two2one), maxDepth);

        tftt::plot::hilbert(formatString("{0}/hilb-t{1}.dat", outDir, two2one));

        for (int l = 0; l < maxDepth; l++) {
            tftt::tikz::meshLayer(formatString("{0}/layer{2}-t{1}.tex", outDir, two2one, l), l);
        }

    }

    tftt::distribute(worldSize);

    tikzCurves(formatString("{0}/hilb-t3-r{1}.tex", outDir, "{0}"));
}
