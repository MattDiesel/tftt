
#include <fstream>
#include <string>
#include <vector>

#include "util/formatstring.h"

#include "../config.h"
#include "../cellref.h"
#include "../tfttops.h"
#include "../iter/all.h"
#include "../fttcore.h"
#include "../structure/tree.h"

#include "serialize.h"


namespace tftt {


void saveTree(std::string fname)
{
    std::ofstream ofs(fname, std::ios::binary);
    saveTree(ofs);
}


template<typename T>
void writeVal(std::ostream& os, T t)
{
    os.write(reinterpret_cast<const char*>(&t), sizeof(T));
}


void writeCell(std::ostream& os, cell_t const& cl)
{
    writeVal(os, cl.id().id);
    writeVal(os, cl.rank());
    writeVal(os, cl.data());
}


void writeHeader(std::ostream& os)
{
    char magic[6] = "LT";
    magic[2] = 3;
    magic[3] = 0;
    magic[4] = DIM;
    magic[5] = sizeof(ident_t);
    os.write(magic, 6);

    uint16_t i = sizeof(data_t);
    writeVal(os, i);

    writeVal<double>(os, 0.0); // Origin fixed at 0.0, 0.0
    writeVal<double>(os, 0.0);

    writeVal(os, gtree.size[0]);
    writeVal(os, gtree.size[1]);
}


void saveTree(std::ostream& os)
{
    os.imbue(std::locale::classic());

    writeHeader(os);

    for (auto const& cl : curve) {
        writeCell(os, cl);
    }
}


void addGhosts(std::set<cell_t>& ghosts, cell_t cl, node_t node)
{
    // Ghosts are any cells required by poisson coefficients

    TreeCell& tc = cl.group->cells[cl.index];

    for (int n = 0; n < tc.poisNgbC; n++) {
        if (tc.poisNgb[n].rank() != node && !tc.poisNgb[n].isBoundary()) {
            ghosts.insert(tc.poisNgb[n]);
        }
    }
}


void splitToDisk(std::string fnameFmt)
{
    int node = 0;
    std::ofstream ofs(::util::formatString(fnameFmt, 0), std::ios::binary);
    writeHeader(ofs);

    std::set<cell_t> ghosts;
    cell_t ngb;

    for (auto const& cl : curve) {
        if (cl.rank() != node) {

            // Write all ghosts
            for (auto& gh : ghosts) {
                writeCell(ofs, gh);
            }
            ghosts.clear();

            ofs.close();
            node++;
            ofs.open(::util::formatString(fnameFmt, node), std::ios::binary);
            writeHeader(ofs);
        }

        writeCell(ofs, cl);

        addGhosts(ghosts, cl, node);
    }

    for (auto& gh : ghosts) {
        writeCell(ofs, gh);
    }
    ghosts.clear();

    ofs.close();
}


template<typename T>
void readVal(std::istream& ist, T& ret)
{
    ist.read(reinterpret_cast<char*>(&ret), sizeof(T));
}


template<typename T>
T& vectorSet(std::vector<T>& vt, unsigned int index)
{
    if (index < vt.size()) {
        return vt[index];
    }

    vt.reserve(index+1);
    while (index+1 >= vt.size()) {
        vt.push_back(T());
    }

    return vt[index];
}


void loadTree(std::string fname, int n)
{
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

    // std::cout << "Dim = " << dim << "\n";
    // std::cout << "IdLen = " << idlen << "\n";
    // std::cout << "DataLen = " << datalen << "\n";

    // Origin not used
    // double x = *reinterpret_cast<double*>(&header[8]);
    // double y = *reinterpret_cast<double*>(&header[16]);
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
        readVal(ifs, toInsert);

        newcl = insert(toInsert);
        if (!newcl.isValid()) {
            throw std::runtime_error("Error loading tree.");
        }

        readVal(ifs, newcl.rank());
        readVal(ifs, newcl.data());
    }

    // Check
    if (n != -1) {
        bool first = false;

        gtree.cactive = 0;

        #ifdef TFTT_DEBUG
        bool nonRank = false;
        #endif

        for (auto& cl : curve) {
            if (cl.rank() == gtree.rank) {
                if (!first) {
                    first = true;
                    gtree.firstActive = cl;
                }
                #ifdef TFTT_DEBUG
                else {
                    if (nonRank) {
                        throw std::runtime_error("Non-contigous range");
                        nonRank = false;
                    }
                }
                #endif

                gtree.cactive++;
                gtree.lastActive = cl;

            }
            else {
                #ifdef TFTT_DEBUG
                if (first) {
                    nonRank = true;
                }
                #endif

                if (cl.rank() != -1) {
                    vectorSet(gtree.ghosts, cl.rank()).insert(cl);
                }
            }
        }
    }

    // Calc border, and fill raw ghost
    for (unsigned int r = 0; r < gtree.ghosts.size(); r++) {
        gtree.rawGhosts.push_back(std::vector<cell_t>());

        for (auto& gh : gtree.ghosts[r]) {
            gtree.rawGhosts[r].push_back(gh);

            // Border cells are any required by the ghosts poisson coefs.
            calcFaceCoefs(gh);

            TreeCell& tc = gh.group->cells[gh.index];

            for (int p = 0; p < tc.poisNgbC; p++) {
                if (tc.poisNgb[p].rank() == n) {
                    vectorSet(gtree.borders, gh.rank()).insert(tc.poisNgb[p]);
                }
            }
        }
    }

    // Fill raw
    for (unsigned int r = 0; r < gtree.borders.size(); r++) {
        gtree.rawBorders.push_back(std::vector<cell_t>());

        for (auto& br : gtree.borders[r]) {
            gtree.rawBorders[r].push_back(br);
        }
    }


    // Fill empty data array
    size_t szg, szb;
    for (unsigned int r = 0; r < gtree.ghosts.size(); r++) {
        szg = gtree.ghosts[r].size();
        szb = gtree.borders[r].size();

        if (szg) {
            gtree.ghostData.push_back(new data_t[szg]);
            gtree.ghostAdaptVectors.push_back(new uint32_t[szg]);
        }
        else {
            gtree.ghostData.push_back(nullptr);
            gtree.ghostAdaptVectors.push_back(nullptr);
        }

        if (szb) {
            gtree.borderData.push_back(new data_t[szb]);
            gtree.borderAdaptVectors.push_back(new uint32_t[szb]);
        }
        else {
            gtree.borderData.push_back(nullptr);
            gtree.borderAdaptVectors.push_back(nullptr);
        }
    }
}


} // namespace tftt
