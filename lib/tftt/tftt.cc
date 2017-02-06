
#include <fstream>
#include <stdexcept>
#include <set>

#include "tftt.h"
#include "tree.h"
#include "treegroup.h"
#include "cellref.h"
#include "gray.h"


namespace tftt {


fnDataNorm dataNorm;


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


cell_t atPos(double pos[DIM]) {
    TreeGroup* gr = gtree.root;

    int n;
    bool loop;
    do {
        n = 0;
        loop = false;
        for (auto& cl : *gr) {
            if (cl.containsPoint(pos)) {
                if (cl.hasChildren()) {
                    gr = cl.children();
                    loop = true;
                    break;
                }
                return CellRef(gr, n);
            }
            n++;
        }
    } while (loop);

    return CellRef();
}

cell_t atVertex(int v) {
    CellRef c(gtree.root, v);

    while (c.hasChildren()) {
        c = c.child(v);
    }

    return c;
}


void refine(CellRef cl) {
    cl.group->cells[cl.index].children = new TreeGroup(cl);
    gtree.ccells += (1 << DIM) - 1;
}


void coarsen(CellRef cl) {

    cell_t nbc;
    for (int nb = 0; nb < 2*DIM; nb++) {
        if (cl.group->neighbours[nb] == cl) {
            cl.group->neighbours[nb] == cl.parent();
        }
    }

    delete cl.group->cells[cl.index].children;
    cl.group->cells[cl.index].children = nullptr;
    gtree.ccells -= (1 << DIM) - 1;
}





void twoToOne_Add(std::set<CellRef>& ls, CellRef cl, CellRef from) {
    int lvl = cl.level();
    CellRef nb;
    for (int n = 0; n < 2*DIM; n++) {
        nb = cl.neighbour(n);
        if (nb == from)
            continue;

        if (nb.isBoundary()) {
            continue;
        }

        if (nb.level() < lvl) {
            ls.insert(nb);
            twoToOne_Add(ls, nb, cl);
        }
    }
}

void twoToOne(CellRef cl) {
    std::set<CellRef> refList;
    twoToOne_Add(refList, cl, CellRef());

    for (auto& cr : refList) {
        refine(cr);
    }
}

std::set<CellRef> adaptList;

void adaptBegin() {
    adaptList.clear();
}
void adaptAdd(CellRef cr) {
    adaptList.insert(cr);
    twoToOne_Add(adaptList, cr, CellRef());
}
bool adaptCommit() {
    for (auto& cr : adaptList) {
        if (!cr.hasChildren())
            refine(cr);
    }
    return adaptList.empty();
}

void adaptAddCoarsen(CellRef cr) {
    cr.group->flaggedForCoarsening = true;
    adaptList.insert(cr);
}
bool adaptCommitCoarsen() {
    for (auto& cr : adaptList) {
        coarsen(cr);
    }
    return adaptList.empty();
}


} // namespace tftt
