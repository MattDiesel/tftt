
#include <iostream>
#include <fstream>
#include <algorithm>

#include "util/formatstring.h"

#include "tftt/tftt.h"
#include "tftt/io/tikz.h"


using namespace util; // for formatString


double epsilonRef = 0.1;
double xshift = -2.0;
double yshift = -1.5;


int escape_iterations(double re, double im, int max_iterations)
{
    double xre = 0.0, xim = 0.0, rtemp, xre2 = 0.0, xim2 = 0.0;
    int i = max_iterations;

    do {
        rtemp = xre2 - xim2 + re;
        xim = 2*xre*xim + im;
        xre = rtemp;
        xre2 = rtemp*rtemp;
        xim2 = xim*xim;
        --i;
    }
    while ( (xre2 + xim2 < 4.0f) && i);

    return i;
}

double mbrot(double re, double im)
{
    return escape_iterations(re, im, 255) / 255.0;
}


void refineLayer()
{
    double eps;

    tftt::adaptSwBegin();

    for (auto& cl : tftt::curve) {
        for (auto& nb : cl.neighbours()) {
            if (nb.isBoundary()) continue;

            eps = cl->cc - nb.avrChildren([](tftt::cell_t const& nb) {
                return nb->cc;
            });
            if (eps < -epsilonRef || eps > epsilonRef) {
                tftt::adaptSwSetRefine(cl);
                break;
            }
        }
    }

    tftt::adaptSwCommit();

    for (auto& cl : tftt::leaves) {
        cl->cc = mbrot(cl.centre(0)+xshift, cl.centre(1)+yshift);
    }
}


int main(int argc, char* argv[])
{
    // Defaults:
    int minDepth = 2, maxDepth = 12;
    tftt::options.two2oneFlag = 3;
    std::string outDir = "analysis";

    constexpr int steps = 12;
    constexpr int numWorlds = 7;

    double xshiftS = -2.5;
    double yshiftS = -2;
    double xshiftStep = 0.1;
    double yshiftStep = 0.1;

    std::vector<int> worldSizes {20,22,24,26,28,30,32};

    double summary[steps][numWorlds];

    int firstStep = 12;
    int ghostSum;
    for (int st = 0; st < steps; st++) {
        std::cout << "Step: " << st << std::endl;
        xshift = xshiftS+xshiftStep*st;
        yshift = yshiftS+yshiftStep*st;

        // Init tree to min depth
        tftt::init(3.0, 3.0);

        std::cout << "Refining to minDepth = " << minDepth << std::endl;
        for (int i = 1; i <= minDepth; i++) {
            for (auto& cl : tftt::leaves) {
                tftt::refine(cl);
            }
        }

        std::cout << "Refining to geometry with maxDepth = " << maxDepth << std::endl;

        for (auto& cl : tftt::leaves) {
            cl->cc = mbrot(cl.centre(0)+xshift, cl.centre(1)+yshift);
        }

        std::ofstream ofs(formatString("{0}/mbrot{1}.dat", outDir, st+firstStep));

        ofs << "Cells ";
        for (auto worldSize : worldSizes) {
            ofs << "W=" << worldSize << " ";
        }
        ofs << "\n";

        for (int d = minDepth; d <= maxDepth; d++) {
            std::cout << "At Depth " << d; std::cout.flush();
            refineLayer();

            std::cout << "." << std::endl;

            for (auto& cl : tftt::curve) {
                tftt::calcFaceCoefs(cl);
            }

            std::cout << "\tCell Count: " << tftt::gtree.ccells << std::endl;

            ofs << tftt::gtree.ccells << " ";

            std::set<tftt::cell_t> ghosts;
            int worldSize;
            for (int w = 0; w < numWorlds; w++) {
                worldSize = worldSizes[w];

                std::cout << "\tWorld Size = " << worldSize; std::cout.flush();
                tftt::distribute(worldSize);

                ghosts.clear();
                ghostSum = 0;

                int node = 0;
                for (auto& cl : tftt::curve) {
                    if (cl.rank() != node) {
                        // std::cout << "\t\tRank " << node << ", ghosts: " << ghosts.size() << "\n";
                        ghostSum += ghosts.size();
                        ghosts.clear();
                        node++;
                    }

                    for (auto& p : cl.poissonNeighbourhood()) {
                        if (p.isBoundary()) continue;

                        if (p.rank() != node) {
                            ghosts.insert(p);
                        }
                    }
                }

                // std::cout << "\t\tRank " << node << ", ghosts: " << ghosts.size() << "\n";
                std::cout << ", Ghosts = " << ghostSum << std::endl;
                summary[st][w] = double(ghostSum) / tftt::gtree.ccells;

                ofs << ghostSum << " ";
            }

            ofs << "\n";
        }

        tftt::reset();
    }

    std::ofstream sm(formatString("{0}/summary2.dat", outDir));
    for (int w = 0; w < numWorlds; w++) {
        sm << worldSizes[w] << " ";
        for (int s = 0; s < steps; s++) {
            sm << summary[s][w] << " ";
        }

        sm << "\n";
    }
}

