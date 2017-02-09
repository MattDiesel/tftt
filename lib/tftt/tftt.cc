
#include <fstream>
#include <stdexcept>
#include <set>
#include <cmath>

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


cell_t findmax(fnData dt, double* maxValRet) {
    cell_t max;
    double maxVal = 0.0;
    double val;
    for (auto& cl : leaves) {
        val = dt(cl.data());
        if (val > maxVal) {
            max = cl;
            maxVal = val;
        }
    }

    if (maxValRet) *maxValRet = maxVal;
    return max;
}


double interpChild(cell_t cl, int ch, int forDir, fnData dt) {
    cell_t forNbCl = cl.neighbour(forDir);

    // if (!forNbCl.hasChildren()) {
    //     return cl.data();
    // }

    forNbCl = forNbCl.child(ch ^ (1 << (forDir >> 1)));

    int awayDir;
    cell_t awayNb;

    double ret = dt(cl.data());
    double intDir;

    for (int d = 0; d < DIM; d++) {
        if (d == forDir >> 1) continue; // Interpolate in desired direction last

        awayDir = (d << 1);
        if (ch & awayDir) awayDir++;

        awayNb = cl.neighbour(awayDir);

        if (awayNb.isBoundary()) {
            continue; // Will just interpolate to the same.
        }

        if (awayNb.hasChildren()) {
            // Interpolate between boundary children
            // Todo: Generalise
            int c1 = ch ^ (1 << d);
            int c2 = c1 ^ (1 << (forDir >> 1));
            intDir = (dt(awayNb.child(c1).data())
                      + dt(awayNb.child(c2).data()))*0.5;

            ret += (intDir - ret) / 3.0;
        }
        else {
            intDir = dt(awayNb.data());
            ret += (intDir - ret) / 4.0;
        }
    }

    // Finally interpolate in desired direction
    intDir = dt(forNbCl.data());
    ret += (intDir - ret) / 3.0;

    return ret;
}

double interpALEVertex(cell_t cl, int v, fnData dt) {
    cell_t ngb = cl.neighbour(v | 1);
    cell_t c2;
    if (ngb.isBoundary()) {
        // Boundary condition?
        throw std::runtime_error("Not Implemented: Automatic BC.");
    }
    else if (ngb.hasChildren()) {
        // Neighbour more refined, use child vertex
        ngb = ngb.child(0);
        return dt(ngb.data());
    }
    else if (ngb.level() < cl.level()) {
        // Neighbour less refined, average two vertices
        // Todo: Generalise to 3d
        c2 = ngb.neighbour(v ^ 2);
        return (dt(ngb.data()) + dt(c2.data()))*0.5;
    }
    else {
        // Same level.
        return dt(ngb.data());
    }
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
        nbc = cl.group->neighbours[nb];
        if (nbc.isBoundary()) continue;

        if (nbc.group->neighbours[nb ^ 1].group == cl.children()) {
            nbc.group->neighbours[nb ^ 1] == cl;
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
