
#include <string>
#include <fstream>
#include <ostream>

#include <iostream> // Todo: Remove

#include "tftt.h"
#include "tree.h"
#include "gray.h"
#include "formatstring.h"

#include "tfttio.h"


namespace tftt {


void drawMesh(std::string fname) {
    std::ofstream ofs(fname);
    drawMesh(ofs);
}

void drawMesh(std::ostream& os) {
    for (auto& c : leaves) {
        drawCell(os, c);
    }
}

void drawCell(std::ostream& os, cell_t const& c) {
    for (int v = 0; v < 1<<DIM; v++) {
        os << c.vertex(utils::toGray(v), 0);
        for (int d = 1; d < DIM; d++) {
            os << " " << c.vertex(utils::toGray(v), d);
        }
        os << "\n";
    }

    os << c.vertex(0, 0);
    for (int d = 1; d < DIM; d++) {
        os << " " << c.vertex(0, d);
    }
    os << "\n\n";
}


void drawPartialMesh(std::string fname, cell_t from, cell_t to) {
    std::ofstream ofs(fname);
    drawPartialMesh(ofs, from, to);
}

void drawPartialMesh(std::string fname) {
    std::ofstream ofs(fname);
    drawPartialMesh(ofs, gtree.firstActive, gtree.lastActive);
}

void drawPartialMesh(std::ostream& os, cell_t from, cell_t to) {
    auto bgn = tagCurve::curve_iterator(from);
    auto end = tagCurve::curve_iterator(to);
    end++;

    for (auto& cl = bgn; cl != end; cl++) {
        drawCell(os, *cl);
    }
}


void drawGhosts(std::string fname) {
    std::ofstream ofs(fname);
    drawGhosts(ofs);
}

void drawGhosts(std::ostream& os) {
    for (auto& c : gtree.ghosts) {
        drawCell(os, c);
    }
}


void drawCurve(std::string fname) {
    std::ofstream ofs(fname);
    drawCurve(ofs);
}

void drawCurve(std::ostream& os) {
    for (auto& c : curve) {
        os << c.centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << c.centre(d);
        }
        
        os << "\n";
    }
}


void drawPartialCurve(std::string fname) {
    std::ofstream ofs(fname);
    drawPartialCurve(ofs);
}

void drawPartialCurve(std::ostream& os) {
    for (auto& c : activecurve) {
        os << c.centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << c.centre(d);
        }
        
        os << "\n";
    }
}


void drawBoundaries(std::string fname) {
    std::ofstream ofs(fname);
    for (int b = 0; b < 2*DIM; b++) {
        drawBoundary(ofs, b);
    }
}

void drawBoundary(std::ostream& os, int b) {
    for (auto& c : boundaryCells(b)) {
        drawCell(os, c);
    }
}


void drawMatrix(std::string fname, int imgW, int imgH, fnDataNorm dataNorm) {
    std::ofstream ofs(fname, std::ios::binary);
    drawMatrix(ofs, imgW, imgH, dataNorm);
}

void drawMatrix(std::ostream& os, int imgW, int imgH, fnDataNorm dataNorm) {
    unsigned char* bmp = new unsigned char[imgW*imgH];

    int x1,y1,w,h,x,y;
    for (auto& cell : leaves) {
        x1 = int((cell.origin(0)) * (imgW-1) / gtree.size[0]);
        y1 = int((cell.origin(1)) * (imgH-1) / gtree.size[1]);
        w = int(cell.size(0) * (imgW-1) / gtree.size[0])+1;
        h = int(cell.size(1) * (imgH-1) / gtree.size[1])+1;

        for (y = y1; y < y1+h; y++) {
            for (x = x1; x < x1 + w; x++) {
                bmp[y*imgW+x] = (unsigned char)(dataNorm(cell.data(), 255));
            }
        }
    }

    os << "P5 " << imgW << " " << imgH << " " << 255 << "\n";

    for (y = 0; y < imgH; y++) {
        for (x = 0; x < imgW; x++) {
            os << (unsigned char)(255-bmp[y*imgW+x]);
        }
    }
}


void saveTree(std::string fname) {
    std::ofstream ofs(fname, std::ios::binary);
    saveTree(ofs);
}

template<typename T>
void writeVal(std::ostream& os, T t) {
    os.write(reinterpret_cast<const char*>(&t), sizeof(T));
}

// void writeDouble(std::ostream& os, double d) {
//     os.write(reinterpret_cast<const char*>(&d), 8);
// }

// void writeIdent(std::ostream& os, ident_t idt) {
//     os.write(reinterpret_cast<const char*>(idt), sizeof(ident_t));
// }

void writeCell(std::ostream& os, cell_t const& cl) {
    // os.write(reinterpret_cast<const char*>(&cl.id().id), sizeof(ident_t));
    // os.write(reinterpret_cast<const char*>(&cell->processor), 1);
    // os.write(reinterpret_cast<const char*>(&cell->data), sizeof(data_t));

    writeVal(os, cl.id().id);
    writeVal(os, cl.rank());
    writeVal(os, cl.data());
}

void writeHeader(std::ostream& os) {
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


void saveTree(std::ostream& os) {
    os.imbue(std::locale::classic());

    writeHeader(os);

    for (auto const& cl : curve) {
        writeCell(os, cl);
    }
}


void addChildren(std::set<cell_t>& ghosts, cell_t ngb, node_t node) {
    for (auto& ngbCh : *ngb.children()) {
        if (ngbCh.hasChildren()) {
            addChildren(ghosts, ngbCh, node);
        }
        else if (ngbCh.rank() != node) {
            if (!ngb.isBoundary())
                ghosts.insert(ngbCh);
        }
    }
}

void splitToDisk(std::string fnameFmt) {
    int node = 0;
    std::ofstream ofs(tftt::utils::formatString(fnameFmt, 0), std::ios::binary);
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
            ofs.open(tftt::utils::formatString(fnameFmt, node), std::ios::binary);
            writeHeader(ofs);
        }

        writeCell(ofs, cl);

        for (int nb = 0; nb < 2*DIM; nb++) {
            ngb = cl.neighbour(nb);
            if (ngb.isBoundary()) continue;

            if (ngb.hasChildren()) {
                addChildren(ghosts, ngb, node);
            }
            else if (ngb.rank() != node) {
                ghosts.insert(ngb);
            }
        }
    }

    for (auto& gh : ghosts) {
        writeCell(ofs, gh);
    }
    ghosts.clear();

    ofs.close();
}


template<typename T>
void readVal(std::istream& ist, T& ret) {
    ist.read(reinterpret_cast<char*>(&ret), sizeof(T));
}


void loadTree(std::string fname, int n) {
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


    if (n != -1) {
        bool first = false;
        bool nonRank = false;
        for (auto& cl : curve) {
            if (cl.rank() == gtree.rank) {
                if (!first) {
                    first = true;
                    gtree.firstActive = cl;
                }
                else {
                    if (nonRank) {
                        std::cout << "Warning: Non-contigous range\n";
                    }
                }

                gtree.lastActive = cl;
            }
            else {
                if (first) {
                    nonRank = true;
                }

                if (cl.rank() != -1) {
                    gtree.ghosts.insert(cl);
                }
            }
        }
    }
}


} // namespace tftt
