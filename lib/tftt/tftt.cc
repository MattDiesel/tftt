
#include <fstream>
#include <stdexcept>

#include "tftt.h"
#include "tree.h"
#include "cellref.h"
#include "gray.h"


namespace tftt {


void init(double w, double h) {
    if (gtree.root)
        throw std::runtime_error("Attempt to initialise tree twice.");

    gtree.ccells = 1<<DIM;

    gtree.size[0] = w;
    gtree.size[1] = h;

    // Init top level cells
    gtree.root = new TreeGroup();
    // TODO: TFTT
    // gtree.root->next = nullptr;
    // gtree.root->prev = nullptr;
}

void reset() {
    delete gtree.root;
    gtree.root = nullptr;
}


cell_t find(ident_t idt) {
    TreeGroup* gr = gtree.root;

    int n = 0;
    while (gr) {
        if (gr->id == idt.group()) {
            return CellRef(gr, (int)idt.orthant());
        }

        gr = gr->cells[idt.orthant(n++)].children;
    }

    return CellRef();
}


void refine(CellRef cl) {
    cl.group->cells[cl.index].children = new TreeGroup(cl);
}


void twoToOne(CellRef cl) {
    int lvl = cl.level();
    for (auto& nb : cl.children()->neighbours) {
        if (nb.isBoundary()) {
            continue;
        }

        if (nb.level() < lvl) {
            refine(nb);
            twoToOne(nb);
        }
    }
}

void drawMesh(std::string fname) {
    std::ofstream ofs(fname);
    drawMesh(ofs);
}

void drawMesh(std::ostream& os) {
    for (auto& c : leaves) {
        for (int v = 0; v < 1<<DIM; v++) {
            os << c.vertex(utils::toGray(v), 0);
            for (int d = 1; d < DIM; d++) {
                os << " " << c.vertex(utils::toGray(v), d);
            }
            os << "\n";
        }
        os << "\n";
    }
}


void drawMatrix(std::string fname) {
    std::ofstream ofs(fname, std::ios::binary);
    drawMatrix(ofs);
}

void drawMatrix(std::ostream& os) {
    int imgW = 1024;
    int imgH = 512;

    unsigned char* bmp = new unsigned char[imgW*imgH];

    int x1,y1,w,h,x,y;
    for (auto& cell : leaves) {
        x1 = int((cell.origin(0)) * (imgW-1) / gtree.size[0]);
        y1 = int((cell.origin(1)) * (imgH-1) / gtree.size[1]);
        w = int(cell.size(0) * (imgW-1) / gtree.size[0])+1;
        h = int(cell.size(1) * (imgH-1) / gtree.size[1])+1;

        for (y = y1; y < y1+h; y++) {
            for (x = x1; x < x1 + w; x++) {
                bmp[y*imgW+x] = (unsigned char)(cell.data() * 255.0);
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


} // namespace tftt
