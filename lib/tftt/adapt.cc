
#include <set>

#include "config.h"
#include "structure/treegroup.h"
#include "structure/treecell.h"
#include "structure/tree.h"
#include "cellref.h"
#include "iter/all.h"
#include "fttcore.h"

#include "adapt.h"


namespace tftt {


void twoToOne_Add(std::set<CellRef, CellRef::less>& ls, CellRef cl, CellRef from)
{
    if (!options.two2oneFlag) return; // No 2-2-1

    int lvl = cl.level();
    CellRef nb;

    if (lvl == 0) return;

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

    std::set<CellRef, CellRef::less> refList;
    twoToOne_Add(refList, cl, CellRef());

    for (auto& cr : refList) {
        if (!cr.hasChildren())
            refine(cr);
    }
}


std::set<CellRef, CellRef::less> adaptList;

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

        if (cl.parent().isValid()) {
            tc = &cl.parent().group->cells[cl.parent().index];
            tc->adaptFlags = AF_NoAction;
        }
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
    adaptList.clear();

    for (auto cl = leaves.begin(); cl != leaves.end(); cl++) {
        tc = &cl->group->cells[cl->index];

        if (tc->adaptFlags == AF_Refine) {
            adaptList.insert(*cl);
        }
        else if (tc->adaptFlags == AF_Coarsen) {
            // Coarsen the parent
            CellRef pr = cl->parent();
            coarsen(pr);
            cl = tagLeaves::leaf_iterator(pr);
        }
    }

    for (auto& cr : adaptList) {
        refine(cr);
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
                // return; // I'm sure this should work.
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
