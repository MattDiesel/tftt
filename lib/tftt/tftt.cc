
#include <fstream>
#include <stdexcept>
#include <set>
#include <cmath>
#include <iostream> // TODO: Remove

#include "util/formatstring.h"

#include "tftt.h"
#include "tree.h"
#include "treegroup.h"
#include "cellref.h"
#include "gray.h"
#include "hilbert.h"


namespace tftt {


TFTTOPTIONS options;


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

    gtree.firstActive = CellRef(gtree.root, hilbChild(0, 0));
    gtree.lastActive = CellRef(gtree.root, hilbChild(0, (1<<DIM)-1));

    cell_t cl;
    cell_t bch;
    TreeGroup* newGrp;
    for (int b = 0; b < 2*DIM; b++) {
        cl = CellRef(gtree.boundGroups, b);

        newGrp = cl.group->cells[cl.index].children = new TreeGroup(cl);

        newGrp->boundary = b;
        newGrp->neighbours[b ^ 1] = CellRef(gtree.root, -1);
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


bool checkAround(cell_t cl, int dist, fnCheckCell check)
{
    cell_t tmp;
    for (int nb = 0; nb < 2*DIM; nb++) {
        tmp = cl;
        for (int p = 0; p < dist; p++) {
            tmp = tmp.neighbour(nb);
            if (tmp.isBoundary()) break;
            if (!check(tmp)) return false;
        }
    }

    return true;
}


bool findAround(cell_t cl, int dist, fnCheckCell check)
{
    cell_t tmp;
    for (int nb = 0; nb < 2*DIM; nb++) {
        tmp = cl;
        for (int p = 0; p < dist; p++) {
            tmp = tmp.neighbour(nb);
            if (tmp.isBoundary()) break;
            if (check(tmp)) return true;
        }
    }

    return false;
}


double interpFace(cell_t cl, int fc, fnData dt)
{
    cell_t ngb = cl.neighbour(fc);
    int clLvl = cl.level();
    int ngbLvl = ngb.level();
    double ret;

    if (ngb.hasChildren()) {
        // Neighbour more refined
        // Average children on that face
        // Todo: Generalise to 3d
        // std::cout << "2\n";

        ret = 0.0;
        for (int n = 0; n < 2; n++) {
            ret += dt(ngb.childOnFace(fc ^ 1, n).data()) / 3.0;
        }

        ret += dt(cl.data()) / 3.0;
    }
    else if (clLvl == ngbLvl) {
        // Same level
        // std::cout << "1\n";

        ret = 0.5*(dt(cl.data()) + dt(ngb.data()));
    }
    else {
        // Neighbour less refined

        int awayDir = (fc ^ 2) & ~1;
        if (cl.index & (1<<(awayDir >> 1))) awayDir ^= 1;

        cell_t ngbngb = ngb.neighbour(awayDir);

        ret = 0.0;
        if (ngbngb.hasChildren()) {
            // std::cout << "4\n";

            for (int n = 0; n < 2; n++) {
                ret += dt(ngbngb.childOnFace(awayDir ^ 1, n).data()) / 18.0;
            }
            ret += dt(ngb.data()) * 2.0 / 9.0;
            ret += dt(cl.data()) * 2.0 / 3.0;
        }
        else {
            // std::cout << "3\n";

            ret += dt(ngbngb.data()) / 12.0;
            ret += dt(ngb.data()) / 4.0;
            ret += dt(cl.data()) * 2.0 / 3.0;
        }
    }

    return ret;
}


double interpChild(cell_t cl, int ch, int forDir, fnData dt)
{
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


double interpALEVertex(cell_t cl, int v, fnData dt)
{
    cell_t ngb = cl.neighbour(v | 1);
    cell_t c2;
    if (ngb.hasChildren()) {
        // Neighbour more refined, use child vertex
        ngb = ngb.child(0);
        return dt(ngb.data());
    }
    else if (ngb.level() < cl.level()) {
        // Neighbour less refined, average two vertices
        // Todo: Generalise to 3d
        c2 = ngb.neighbour(v ^ 2);
        if (c2.hasChildren())
            c2 = c2.child(0);

        return (dt(ngb.data()) + dt(c2.data()))*0.5;
    }
    else {
        // Same level.
        return dt(ngb.data());
    }
}


void distribute(int n)
{
    int cpernode = (gtree.ccells / n)+1;

    int cc = cpernode;
    int node = 0;
    for (auto& cl : curve) {
        if (!--cc) {
            cc = cpernode;
            node++;
        }

        cl.rank() = node;
    }
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


void refine(CellRef cl)
{
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
    cell_t nbc;
    for (int nb = 0; nb < 2*DIM; nb++) {
        nbc = cl.group->neighbours[nb];
        if (nbc.isBoundary() && nbc.level() >= cl.level()) {
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


void twoToOne_Add(std::set<CellRef, crless>& ls, CellRef cl, CellRef from)
{
    if (!options.two2oneFlag) return; // No 2-2-1

    int lvl = cl.level();
    CellRef nb;

    if (options.two2oneFlag == -1) {
        for (int n = 0; n < 2*DIM; n++) {
            nb = cl.neighbour(n);
            if (nb == from)
                break;

            if (nb.isBoundary()) {
                break;
            }

            if (nb.level() < lvl) {
                ls.insert(nb);
                twoToOne_Add(ls, nb, cl);
            }
        }

        // Include corners
        cell_t nb2;
        for (int n = 0; n < DIM-1; n++) {
            for (int b = 0; b <= 1; b++) {
                nb = cl.neighbour(n << 1 | b);
                if (nb == from)
                    continue;

                if (nb.isBoundary()) {
                    continue;
                }

                for (int a = 0; a < DIM; a++) {
                    if (a == n) continue;

                    for (int dir = 0; dir <= 1; dir++) {
                        if (nb.level() < lvl
                                && (dir != ((cl.index >> a) & 1))) {
                            // If neighbour is at a lower level, its neighbour
                            // won't be on the corner unless it's in the same direction
                            continue;
                        }

                        nb2 = nb.neighbour(a << 1 | dir);

                        if (nb2.level() < lvl) {
                            ls.insert(nb2);
                            twoToOne_Add(ls, nb2, cl);
                        }
                    }
                }
            }
        }
    }
    else {
        // Propagate in each direction p times

        for (int n = 0; n < 2*DIM; n++) {
            nb = cl;

            // Propagate in given direction
            for (int p = 0; p < options.two2oneFlag; p++) {
                nb = nb.neighbour(n);
                if (nb == from)
                    break;

                if (nb.isBoundary()) {
                    break;
                }

                if (nb.level() < lvl) {
                    p++;
                    ls.insert(nb);
                    twoToOne_Add(ls, nb, cl);
                }
            }
        }
    }
}

void twoToOne(CellRef cl)
{
    if (!options.two2oneFlag) return; // No 2-2-1

    std::set<CellRef, crless> refList;
    twoToOne_Add(refList, cl, CellRef());

    for (auto& cr : refList) {
        if (!cr.hasChildren())
            refine(cr);
    }
}


std::set<CellRef, crless> adaptList;

void adaptBegin()
{
    adaptList.clear();
}
void adaptAdd(CellRef cr)
{
    adaptList.insert(cr);
    twoToOne_Add(adaptList, cr, CellRef());
}
bool adaptCommit()
{
    for (auto& cr : adaptList) {
        if (!cr.hasChildren())
            refine(cr);
    }
    return adaptList.empty();
}

void adaptAddCoarsen(CellRef cr)
{
    adaptList.insert(cr);
}
bool adaptCommitCoarsen()
{

    // int inter = 0;
    // Do it reverse order so that 2-1 is maintained.
    for (auto cr = adaptList.rbegin(); cr != adaptList.rend(); cr++) {
        coarsen(*cr);

        // drawCurve(::util::formatString("coarsen.{0}.dat", inter++));
    }
    return adaptList.empty();
}


void adaptSwBegin()
{
    TreeCell* tc;
    for (auto cl : leaves) {
        tc = &cl.group->cells[cl.index];
        tc->adaptFlags = AF_NoAction;

        tc = &cl.parent().group->cells[cl.parent().index];
        tc->adaptFlags = AF_NoAction;
    }

    for (auto& vec : gtree.ghostAdaptVectors) {
        vec = 0;
    }
}


void adaptSwSetFlags(CellRef cl, ADAPTFLAGS af)
{
    TreeCell& tc = cl.group->cells[cl.index];

    if (cl.hasChildren()) {
        if (af == AF_Refine) {
            af = AF_HoldRefined;
            // throw std::logic_error("Cell already refined!");
        }

        for (int ch = 0; ch < 1<<DIM; ch++) {
            cell_t clch = cl.child(ch);
            if (af == AF_HoldRefined && clch.hasChildren()) continue;

            adaptSwSetFlags(cl.child(ch), af);
        }
    }
    else {
        if ((af == AF_Coarsen || af == AF_HoldCoarsened)) {
            if (tc.adaptFlags != AF_NoAction)
                return; // Ignore
        }

        if (tc.adaptFlags == AF_Refine && af == AF_HoldRefined)
            return; // Nothing to do here

        if (tc.adaptFlags == af) return;

        if (tc.adaptFlags != AF_NoAction && tc.adaptFlags != af) {
            if (tc.adaptFlags == AF_HoldRefined && af == AF_Refine) {
            }
            else
                throw std::runtime_error("Conflicting refinement instructions");
        }
    }

    tc.adaptFlags = af;
}

ADAPTFLAGS adaptSwGetFlags(CellRef cl)
{
    TreeCell& tc = cl.group->cells[cl.index];
    return tc.adaptFlags;
}


void adaptSwCommit()
{
    TreeCell* tc;

    for (auto cl = leaves.begin(); cl != leaves.end(); cl++) {
    nextLeaf:
        tc = &cl->group->cells[cl->index];

        if (tc->adaptFlags == AF_Refine) {
            // Store cell and increment iterator before refinement
            CellRef old = *cl;
            cl++;

            refine(old);

            // Skip child processing:
            goto nextLeaf;
        }
        else if (tc->adaptFlags == AF_Coarsen) {
            // Coarsen the parent
            CellRef pr = cl->parent();
            coarsen(pr);
            cl = tagLeaves::leaf_iterator(pr);
        }


    }
}


void adaptSwSetCoarsen(CellRef cl)
{
    if (cl.hasChildren()) {
        for (int ch = 0; ch < 1<<DIM; ch++) {
            switch (adaptSwGetFlags(cl.child(ch))) {
                case AF_Refine:
                case AF_HoldRefined:
                    return;
                default:
                    break;
            }
        }
        adaptSwSetFlags(cl, AF_Coarsen);
    }
    else {
        adaptSwSetFlags(cl, AF_HoldCoarsened);
    }
}


void adaptSwPropogateLevel(CellRef cl, int dir, int lvl)
{
    for (int d = lvl; d > 1; d--) {
        for (int p = 0; p < options.two2oneFlag; p++) {
            cl = cl.neighbour(dir);
            if (cl.isBoundary()) return;

            if (cl.hasChildren()) {
                // Do nothing
            }
            else if (cl.level() == d) {
                adaptSwSetFlags(cl.parent(), AF_HoldRefined);
                adaptSwPropogateLevel(cl, dir ^ 2, d-1);
                adaptSwPropogateLevel(cl, dir ^ 3, d-1);
            }
            else if (cl.level() < d) {
                p++;
                adaptSwSetFlags(cl, AF_Refine);
                adaptSwPropogateLevel(cl, dir ^ 2, d-1);
                adaptSwPropogateLevel(cl, dir ^ 3, d-1);
            }
        }
        if (cl.level() == d) cl = cl.parent();
    }
}


void adaptSwSetRefine(CellRef cl)
{
    if (cl.hasChildren()) {
        throw std::logic_error("Cell already refined.");
    }

    adaptSwSetFlags(cl, AF_Refine);
    for (int d = 0; d < 4; d++) {
        adaptSwPropogateLevel(cl, d, cl.level());
    }
}


void adaptSwSetHoldRefined(CellRef cl)
{
    adaptSwSetFlags(cl.parent(), AF_HoldRefined);
    for (int d = 0; d < 4; d++) {
        adaptSwPropogateLevel(cl.parent(), d, cl.level()-1);
    }
}


} // namespace tftt
