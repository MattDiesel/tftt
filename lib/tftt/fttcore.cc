
#include <stdexcept>

#include "config.h"
#include "hilbert.h"
#include "cellref.h"
#include "tree.h"
#include "treegroup.h"

#include "fttcore.h"


namespace tftt {


void init(double w, double h)
{
    if (gtree.root)
        throw std::runtime_error("Attempt to initialise tree twice.");

    gtree.rank = -1;

    gtree.ccells = 1<<DIM;

    gtree.size[0] = w;
    gtree.size[1] = h;

    for (int b = 0; b < 2*DIM; b++) {
        gtree.boundGroups = new TreeGroup();
    }

    // Init top level cells
    gtree.root = new TreeGroup();
    cell_t cl;
    cell_t bch;
    TreeGroup* newGrp;


    for (int b = 0; b < 2*DIM; b++) {
        cl = CellRef(gtree.boundGroups, b);
        cl.group->cells[cl.index].children = new TreeGroup(cl);
    }

    for (int b = 0; b < 2*DIM; b++) {
        cl = CellRef(gtree.boundGroups, b);
        newGrp = cl.group->cells[cl.index].children;

        newGrp->boundary = b;
        for (int b2 = 0; b2 < 2*DIM; b2++) {
            if (b2 == (b ^ 1))
                newGrp->neighbours[b2] = CellRef(gtree.root, -1);
            else
                newGrp->neighbours[b2] = CellRef();
        }
        gtree.root->neighbours[b] = CellRef(gtree.boundGroups, b);

        newGrp->id = 0;
        for (int d = 0; d < DIM; d++) {
            newGrp->origin[d] = 0.0;

            if (d*2 == b)
                newGrp->origin[d] -= gtree.size[d];
            else if (d*2 == (b^1))
                newGrp->origin[d] += gtree.size[d];

            // std::cout << "(" << b << "," << d << ") @ " << newGrp->origin[d] << std::endl;
        }
    }

    gtree.first = gtree.firstActive = CellRef(gtree.root, hilbChild(0, 0));
    gtree.last = gtree.lastActive = CellRef(gtree.root, hilbChild(0, (1<<DIM)-1));


    gtree.destroying = false;
}


void reset()
{
    gtree.destroying = true;
    delete gtree.root;
    delete gtree.boundGroups;
    gtree.root = nullptr;
    gtree.ghosts.clear();
}


void refine(CellRef cl)
{
    #ifdef TFTT_DEBUG
    if (cl.isBoundary())
        throw std::runtime_error("Calling refine() on boundary.");
    #endif

    cl.group->cells[cl.index].children = new TreeGroup(cl);
    gtree.ccells += (1 << DIM) - 1;

    // Refine boundary as needed.
    cell_t ngb;
    for (int nb = 0; nb < 2*DIM; nb++) {
        ngb = cl.neighbour(nb);
        if (ngb.isBoundary()) {
            ngb.group->cells[ngb.index].children = new TreeGroup(ngb);
        }
    }
}


void coarsen(CellRef cl)
{
    #ifdef TFTT_DEBUG
    for (int ch = 0; ch < (1<<DIM); ch++) {
        if (cl.child(ch).hasChildren()) {
            throw std::runtime_error("Cell has grandchildren.");
        }
    }
    #endif

    cell_t nbc;
    for (int nb = 0; nb < 2*DIM; nb++) {
        nbc = cl.children()->neighbours[nb];

        if (nbc.isBoundary()) {
            #ifdef TFTT_DEBUG
            if (nbc.level() != cl.level())
                throw std::runtime_error("Boundary level mismatch");
            #endif

            delete nbc.group->cells[nbc.index].children;
            nbc.group->cells[nbc.index].children = nullptr;
        };

        if (nbc.group->neighbours[nb ^ 1].group == cl.children()) {
            nbc.group->neighbours[nb ^ 1] == cl;
        }
    }

    delete cl.group->cells[cl.index].children;
    cl.group->cells[cl.index].children = nullptr;
    gtree.ccells -= (1 << DIM) - 1;
}


cell_t atPos(double pos[DIM])
{
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
    }
    while (loop);

    return CellRef();
}


cell_t atVertex(int v)
{
    CellRef c(gtree.root, v);

    while (c.hasChildren()) {
        c = c.child(v);
    }

    return c;
}


} // namespace tftt
