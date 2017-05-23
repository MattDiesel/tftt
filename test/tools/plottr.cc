

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <map>

#include "tftt/tftt.h"
#include "tftt/structure/tree.h"

#include "util/formatstring.h"

using namespace tftt;

template<typename T>
void readVal2(std::istream& ist, T& ret)
{
    ist.read(reinterpret_cast<char*>(&ret), sizeof(T));
}


void printTree(std::string fname)
{
    std::string tempPath = "./plottrtmp";
    std::map<int, std::ofstream> ofss;

    std::ifstream ifs(fname, std::ios::binary);

    TrHeader trh;
    readVal2(ifs, trh);
    trh.check();

    init(trh.domainSize[0], trh.domainSize[1]);
    gtree.rank = trh.worldRank;

    ident_t toInsert;
    cell_t newcl;

    while (!ifs.eof()) {
        readVal2(ifs, toInsert);

        newcl = insert(toInsert);
        if (!newcl.isValid()) {
            throw std::runtime_error("Error loading tree.");
        }

        std::cout << newcl << " inserted. ";

        readVal2(ifs, newcl.rank());
        std::cout << "Rank=" << int(newcl.rank()) << "\n";

        readVal2(ifs, newcl.data());

        auto search = ofss.find(newcl.rank());

        if (search == ofss.end()) {
            ofss[newcl.rank()] = std::ofstream(::util::formatString("{0}/r{1}.dat", tempPath, (int)newcl.rank()));
            search = ofss.find(newcl.rank());
        }

        drawCell(search->second, newcl);
    }
}


int main(int argc, char const* argv[])
{

    if (argc < 2) {
        std::cout << "Usage: plottr [file]" << std::endl;
        return 1;
    }

    tftt::options.two2oneFlag = 0;

    printTree(argv[1]);

    return 0;
}
