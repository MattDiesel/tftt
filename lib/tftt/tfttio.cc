
#include <string>
#include <fstream>
#include <ostream>

#include <iostream> // Todo: Remove

#include "util/formatstring.h"

#include "tftt.h"
#include "tree.h"
#include "gray.h"

#include "tfttio.h"


namespace tftt {


void drawMeshSub2d(std::ostream& os, TreeGroup* gr, double w, double h, double x, double y)
{
    // if (!gr) {
    os << x << ' ' << y << '\n'
       << (x+w) << ' ' << y << '\n'
       << (x+w) << ' ' << (y+h) << '\n'
       << x << ' ' << (y+h) << '\n'
       << x << ' ' << y << "\n\n";
    // }
    // else {
    if (gr) {
        drawMeshSub2d(os, gr->cells[0].children, w*0.5, h*0.5, x, y);
        drawMeshSub2d(os, gr->cells[1].children, w*0.5, h*0.5, x+w*0.5, y);
        drawMeshSub2d(os, gr->cells[2].children, w*0.5, h*0.5, x, y+h*0.5);
        drawMeshSub2d(os, gr->cells[3].children, w*0.5, h*0.5, x+w*0.5, y+h*0.5);
    }
}


void drawMesh(std::string fname)
{
    std::ofstream ofs(fname);
    drawMesh(ofs);
}

void drawMesh(std::ostream& os)
{
    if (DIM == 2) {
        drawMeshSub2d(os, gtree.root, gtree.size[0], gtree.size[1], 0.0, 0.0);
    }
    else {
        for (auto& c : leaves) {
            drawCell(os, c);
        }
    }
}


void drawCell(std::ostream& os, cell_t const& c)
{
    double vtc[1<<DIM][DIM];

    c.vertices(vtc);

    int vg;
    for (int v = 0; v < 1<<DIM; v++) {
        vg = utils::toGray(v);

        // os << c.vertex(utils::toGray(v), 0);
        os << vtc[vg][0];
        for (int d = 1; d < DIM; d++) {
            // os << " " << c.vertex(vg, d);
            os << ' ' << vtc[vg][d];
        }
        os << '\n';
    }

    // os << c.vertex(0, 0);
    os << vtc[0][0];
    for (int d = 1; d < DIM; d++) {
        // os << " " << c.vertex(0, d);
        os << ' ' << vtc[0][d];
    }
    os << "\n\n";
}


void drawPartialMesh(std::string fname, cell_t from, cell_t to)
{
    std::ofstream ofs(fname);
    drawPartialMesh(ofs, from, to);
}

void drawPartialMesh(std::string fname)
{
    std::ofstream ofs(fname);
    drawPartialMesh(ofs, gtree.firstActive, gtree.lastActive);
}

void drawPartialMesh(std::ostream& os, cell_t from, cell_t to)
{
    auto bgn = tagCurve::curve_iterator(from);
    auto end = tagCurve::curve_iterator(to);
    end++;

    for (auto& cl = bgn; cl != end; cl++) {
        drawCell(os, *cl);
    }
}


void drawGhosts(std::string fname)
{
    std::ofstream ofs(fname);
    drawGhosts(ofs);
}

void drawGhosts(std::ostream& os)
{
    for (auto& c : gtree.ghosts) {
        drawCell(os, c);
    }
}


void drawCurve(std::string fname)
{
    std::ofstream ofs(fname);
    drawCurve(ofs);
}

void drawCurve(std::ostream& os)
{
    for (auto& c : curve) {
        os << c.centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << c.centre(d);
        }

        os << "\n";
    }
}


void drawPartialCurve(std::string fname)
{
    std::ofstream ofs(fname);
    drawPartialCurve(ofs);
}

void drawPartialCurve(std::ostream& os)
{
    for (auto& c : activecurve) {
        os << c.centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << c.centre(d);
        }

        os << "\n";
    }
}


void drawBoundaries(std::string fname)
{
    std::ofstream ofs(fname);
    for (int b = 0; b < 2*DIM; b++) {
        drawBoundary(ofs, b);
    }
}

void drawBoundary(std::ostream& os, int b)
{
    for (auto& c : boundaryCells(b)) {
        drawCell(os, c);
    }
}


void drawMatrix(std::string fname, int imgW, int imgH, fnDataNorm dataNorm)
{
    std::ofstream ofs(fname, std::ios::binary);
    drawMatrix(ofs, imgW, imgH, dataNorm);
}

void drawMatrix(std::ostream& os, int imgW, int imgH, fnDataNorm dataNorm)
{
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


void plotMatrix(std::string fname, fnData dt)
{
    std::ofstream ofs(fname);
    plotMatrix(ofs, dt);
}

void plotMatrix(std::ostream& os, fnData dt)
{
    for (auto& cl : leaves) {
        os << cl.centre(0) << " " << cl.centre(1) << " " << dt(cl.data()) << "\n";
    }
}


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

// void writeDouble(std::ostream& os, double d) {
//     os.write(reinterpret_cast<const char*>(&d), 8);
// }

// void writeIdent(std::ostream& os, ident_t idt) {
//     os.write(reinterpret_cast<const char*>(idt), sizeof(ident_t));
// }

void writeCell(std::ostream& os, cell_t const& cl)
{
    // os.write(reinterpret_cast<const char*>(&cl.id().id), sizeof(ident_t));
    // os.write(reinterpret_cast<const char*>(&cell->processor), 1);
    // os.write(reinterpret_cast<const char*>(&cell->data), sizeof(data_t));

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


void addChildren(std::set<cell_t>& ghosts, cell_t cl, cell_t ngb, int nb, node_t node)
{
    for (auto& ngbCh : *ngb.children()) {
        if (options.ghostsFlag == 0 && ngbCh.neighbour(nb ^ 1) != cl)
            continue; // Minimal - only take bordering children

        if (ngbCh.hasChildren()) {
            addChildren(ghosts, cl, ngbCh, nb, node);
        }
        else if (ngbCh.rank() != node) {
            ghosts.insert(ngbCh);
        }
    }
}


void addGhosts(std::set<cell_t>& ghosts, cell_t cl, node_t node)
{
    cell_t ngb;
    for (int nb = 0; nb < 2*DIM; nb++) {
        ngb = cl.neighbour(nb);
        if (ngb.isBoundary()) continue;

        if (ngb.hasChildren()) {
            addChildren(ghosts, cl, ngb, nb, node);
        }
        else if (ngb.rank() != node) {
            ghosts.insert(ngb);
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

        for (auto& cl : curve) {
            if (cl.rank() == gtree.rank) {
                if (!first) {
                    first = true;
                    gtree.firstActive = cl;
                }
                #ifdef TFTT_DEBUG
                else {
                    if (nonRank) {
                        std::cout << "Warning: Non-contigous range " << cl << "\n";
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
                    gtree.ghosts.insert(cl);
                }
            }
        }
    }
}





void drawPoissonNeighbourhood(std::string fname, cell_t cl)
{
    std::ofstream ofs(fname);
    drawPoissonNeighbourhood(ofs, cl);
}

void drawPoissonNeighbourhood(std::ostream& os, cell_t cl)
{

    TreeCell* tc = &cl.group->cells[cl.index];

    for (int n = 0; n < tc->poisNgbC; n++) {
        drawCell(os, tc->poisNgb[n]);
    }

}





} // namespace tftt
