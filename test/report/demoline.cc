
#include <iostream>
#include <fstream>
#include <algorithm>

#include "util/formatstring.h"

#include "tftt/tftt.h"


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

void tikzCell(std::ostream& os, tftt::cell_t cl)
{
    // os << "\\draw (" << cl.origin(0) << "," << cl.origin(1)
    //    << ") rectangle (" << (cl.origin(0)+cl.size(0)) << "," << (cl.origin(1)+cl.size(1)) << ");" << "\n";
    os << "\\node [draw,shape=rectangle,anchor=center] at ("
       << cl.centre(0) << "," << cl.centre(1) << ") {};\n";
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


void tikzLayer(std::ostream& os, tftt::cell_t cl, int layer)
{
    if (cl.level() == layer) {
        tikzCell(os, cl);
    }
    else if (cl.hasChildren()) {
        for (auto ch : *cl.children()) {
            tikzLayer(os, ch, layer);
        }
    }
}

void tikzLayer(std::ostream& ofs, int layer)
{
    tftt::cell_t root = tftt::cell_t(-1);

    double lw = 0.1;
    double w = tftt::gtree.size[0] / (2 << layer) - 2*lw*0.001;
    double h = tftt::gtree.size[1] / (2 << layer) - 2*lw*0.001;

    ofs << "\\begin{scope}[minimum width=" << w << "cm,minimum height=" << h << "cm,inner sep=0, line width=" << lw << "mm]\n";

    tikzLayer(ofs, root, layer);

    ofs << "\\end{scope}\n";
}

void tikzLayer(std::string fname, int layer)
{
    std::ofstream ofs(fname);
    tikzLayer(ofs, layer);
}

void tikzMesh(std::string fname, int maxDepth)
{
    std::ofstream ofs(fname);

    for (int d = 0; d < maxDepth; d++) {
        tikzLayer(ofs, d);
    }
}

int main(int argc, char* argv[])
{
    // Defaults:
    int minDepth = 0, maxDepth = 5;
    double p[2];

    std::string outDir = "../docs/report/method";


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

        tikzMesh(formatString("{0}/mesh-t{1}.tex", outDir, two2one), maxDepth);

        for (int l = 0; l < maxDepth; l++) {
            tikzLayer(formatString("{0}/layer{2}-t{1}.tex", outDir, two2one, l), l);
        }

        tftt::reset();
    }
}
