
#include <iostream>
#include <fstream>
#include <algorithm>

#include "util/formatstring.h"

#include "tftt/tftt.h"
#include "tftt/io/tikz.h"


using namespace util; // for formatString


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
    int minDepth = 0, maxDepth = 3;
    int worldSize = 3;
    tftt::options.two2oneFlag = 1;

    double p[2] = {0.25, 2.75};

    std::string outDir = "../docs/report/method/gen";



    // Init tree to min depth
    tftt::init(4.0, 4.0);

    tftt::tikz::hilbert(formatString("{0}/tr-hilb-l0.tex", outDir));
    tftt::tikz::morton(formatString("{0}/tr-mort-l0.tex", outDir));

    p[1] = 1;
    tftt::refine(tftt::atPos(p));
    p[1] = 2.75;
    tftt::refine(tftt::atPos(p));

    tftt::tikz::hilbert(formatString("{0}/tr-hilb-l1.tex", outDir));
    tftt::tikz::morton(formatString("{0}/tr-mort-l1.tex", outDir));

    tftt::refine(tftt::atPos(p));

    tftt::tikz::hilbert(formatString("{0}/tr-hilb-l2.tex", outDir));
    tftt::tikz::morton(formatString("{0}/tr-mort-l2.tex", outDir));


    tftt::tikz::hilbert(formatString("{0}/tr-hilb.tex", outDir));
    tftt::tikz::morton(formatString("{0}/tr-mort.tex", outDir));

    tftt::tikz::mesh(formatString("{0}/tr-mesh.tex", outDir), maxDepth);
    // for (int l = 0; l <= maxDepth; l++) {
    //     tftt::tikz::meshLayer(formatString("{0}/tr-layer{1}.tex", outDir, l), l);
    // }
}
