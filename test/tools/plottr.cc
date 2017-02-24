

#include <iostream>
#include <fstream>
#include <stdexcept>

#include "tftt/tftt.h"
#include "tftt/tree.h"

using namespace tftt;

template<typename T>
void readVal2(std::istream& ist, T& ret) {
    ist.read(reinterpret_cast<char*>(&ret), sizeof(T));
}


void printTree(std::string fname, int n) {
    std::ifstream ifs(fname, std::ios::binary);

    // Read and check header
    char header[40];
    ifs.read(reinterpret_cast<char*>(header), 40);

    if (header[0] != 'L' || header[1] != 'T') {
        throw std::invalid_argument(fname);
        throw std::invalid_argument("Input file is not a valid ltree file.");
    }

    int dim = int(header[4]);
    int idlen = int(header[5]);
    uint16_t datalen = *reinterpret_cast<uint16_t*>(&header[6]);

    double w = *reinterpret_cast<double*>(&header[24]);
    double h = *reinterpret_cast<double*>(&header[32]);

    if (header[2] != 3)
        throw std::runtime_error("Major version of file format incompatible");
    if (dim != DIM)
        throw std::runtime_error("Dimension of file tree does not match structure.");
    if (datalen != sizeof(data_t))
        throw std::runtime_error("Data size in file does not match structure.");
    if (idlen != sizeof(ident_t))
        throw std::runtime_error("ID size in file does not match structure.");

    init(w, h);
    gtree.rank = n;

    ident_t toInsert;
    cell_t newcl;

    while (!ifs.eof()) {
        readVal2(ifs, toInsert);

        newcl = insert(toInsert);
        if (!newcl.isValid()) {
            throw std::runtime_error("Error loading tree.");
        }

        readVal2(ifs, newcl.rank());
        readVal2(ifs, newcl.data());

        std::cout << newcl.centre(0);
        for (int d = 1; d < DIM; d++) {
            std::cout << " " << newcl.centre(d);
        }
        
        std::cout << "\n";
    }
}


int main(int argc, char const *argv[])
{

    if (argc < 2) {
        std::cout << "Usage: plottr [file] [n]" << std::endl;
        return 1;
    }

    int n = -1;
    if (argc > 2) {
        n = argv[2][0] - '0';
    }

    tftt::options.two2oneFlag = 0;

    printTree(argv[1], n);

    return 0;
}
