
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


void _drawMeshSub(std::ostream& os, TreeGroup* gr) {
    for (auto& c : *gr) {
        if (c.hasChildren()) 
            _drawMeshSub(os, c.children());
        else {
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
}

void drawMesh(std::string fname) {
    std::ofstream ofs(fname);
    drawMesh(ofs);
}

void drawMesh(std::ostream& os) {
    _drawMeshSub(os, gtree.root);
}



} // namespace tftt
