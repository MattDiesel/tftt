
#include <fstream>
#include <string>
#include <vector>
#include <set>

#include "util/formatstring.h"

#include "../config.h"
#include "../cellref.h"
#include "../tfttops.h"
#include "../iter/all.h"
#include "../fttcore.h"
#include "../structure/tree.h"

#include "serialize.h"


namespace tftt {


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


template<typename T>
void readVal(std::istream& ist, T& ret)
{
    ist.read(reinterpret_cast<char*>(&ret), sizeof(T));
}


TrHeader::TrHeader()
{
    magic = TrHeader::Magic;
    versionMajor = TrHeader::VersionMajor;
    versionMinor = TrHeader::VersionMinor;

    dimensions = DIM;
    identSize = sizeof(ident_t);
    dataSize = sizeof(data_t);

    for (int d = 0; d < DIM; d++) {
        domainOrigin[d] = 0.0;
        domainSize[d] = gtree.size[d];
    }
}


TrHeader::TrHeader(int rank, int world)
    : TrHeader()
{
    worldSize = world;
    worldRank = rank;
}


void TrHeader::check()
{
    if (magic != TrHeader::Magic)
        throw std::runtime_error("Magic Number Mismatch.");
    if (versionMajor != TrHeader::VersionMajor)
        throw std::runtime_error("Major version of file format incompatible");
    if (dimensions != DIM)
        throw std::runtime_error("Dimension of file tree does not match structure.");
    if (identSize != sizeof(ident_t))
        throw std::runtime_error("ID size in file does not match structure.");
    if (dataSize != sizeof(data_t))
        throw std::runtime_error("Data size in file does not match structure.");
}


void saveTree(std::string fname)
{
    std::ofstream ofs(fname, std::ios::binary);
    saveTree(ofs);
}


void saveTree(std::ostream& os)
{
    os.imbue(std::locale::classic());

    TrHeader trh;
    writeVal(os, trh);

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


void saveParTree(std::string fnameFmt, int world)
{
    int node = 0;
    TrHeader trh(0, world);

    std::ofstream ofs(::util::formatString(fnameFmt, 0), std::ios::binary);
    writeVal(ofs, trh);

    std::set<cell_t> ghosts;
    cell_t ngb;

    for (auto const& cl : curve) {
        if (cl.rank() != node) {

            // Write all ghosts
            for (auto& gh : ghosts) {
                writeCell(ofs, gh);
            }
            ghosts.clear();

            // Step onto next file
            ofs.close();

            node++;
            trh.worldRank++;

            ofs.open(::util::formatString(fnameFmt, node), std::ios::binary);
            writeVal(ofs, trh);
        }

        writeCell(ofs, cl);

        addGhosts(ghosts, cl, node);
    }

    // Write remaining ghosts
    for (auto& gh : ghosts) {
        writeCell(ofs, gh);
    }
    ghosts.clear();

    ofs.close();
}


// Deprecated.
// template<typename T>
// T& vectorSet(std::vector<T>& vt, unsigned int index)
// {
//     if (index < vt.size()) {
//         return vt[index];
//     }

//     vt.reserve(index+1);
//     while (index+1 >= vt.size()) {
//         vt.push_back(T());
//     }

//     return vt[index];
// }


void loadParTree(std::string fname)
{
    std::ifstream ifs(fname, std::ios::binary);

    // Check the file header
    TrHeader trh;
    readVal(ifs, trh);
    trh.check();

    init(trh.domainSize[0], trh.domainSize[1]);
    gtree.rank = trh.worldRank;

    // Initialise ghosts/border sets.
    gtree.ghosts.reserve(trh.worldSize);
    gtree.borders.reserve(trh.worldSize);
    for (int r = 0; r < trh.worldSize; r++) {
        gtree.ghosts.push_back(decltype(gtree.ghosts)::value_type());
        gtree.borders.push_back(decltype(gtree.borders)::value_type());
    }

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
    // if (n != -1) {
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
                gtree.ghosts[cl.rank()].insert(cl);
            }
        }
    }
    // }

    // Calc border, and fill raw ghost
    for (unsigned int r = 0; r < gtree.ghosts.size(); r++) {
        gtree.rawGhosts.push_back(std::vector<cell_t>());

        for (auto& gh : gtree.ghosts[r]) {
            gtree.rawGhosts[r].push_back(gh);

            // Border cells are any required by the ghosts poisson coefs.
            calcFaceCoefs(gh);

            TreeCell& tc = gh.group->cells[gh.index];

            for (int p = 0; p < tc.poisNgbC; p++) {
                if (tc.poisNgb[p].rank() == gtree.rank) {
                    gtree.borders[gh.rank()].insert(tc.poisNgb[p]);
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
