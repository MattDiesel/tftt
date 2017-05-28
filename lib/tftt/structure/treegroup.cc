
#include <stdexcept>
#include <iostream>

#include "../config.h"
#include "../cellref.h"
#include "tree.h"
#include "../hilbert.h"
#include "treeface.h"
#include "treevertex.h"
#include "groupmemory.h"

#include "treegroup.h"


namespace tftt {

TreeGroup::TreeGroup()
    : parent()
{
    for (int i = 0; i < (1<<DIM); i++) {
        cells[i].children = nullptr;
        cells[i].rank = -1;
        cells[i].index = i;
    }

    id = 0;
    for (int i = 0; i < DIM; i++) {
        origin[i] = 0.0;
    }

    for (int i = 0; i < DIM*2; i++) {
        neighbours[i] = CellRef(gtree.boundGroups, i);
    }

    // Thread
    orientation = 0;
    next = CellRef();
    prev = CellRef();

    #ifdef TFTT_FACES
    // Boundaries
    // Internal
    constexpr int pairs[4][4] = {
        {0, 1, 1, 0},
        {0, 2, 2, 3},
        {1, 2, 3, 3},
        {2, 1, 3, 0}
    };

    for (int i = 0; i < 4; i++) {
        cells[pairs[i][0]].faces[pairs[i][1]] =
            cells[pairs[i][2]].faces[pairs[i][3]] =
                new TreeFace(CellRef(this, pairs[i][0]), CellRef(this, pairs[i][2]), pairs[i][1] >> 1);
    }

    // Add external faces
    constexpr int facecells[4][2] = {
        {0, 2},
        {1, 3},
        {0, 1},
        {2, 3}
    };

    cell_t cl;
    for (int n = 0; n < 2*DIM; n++) {
        // Add faces
        for (int i = 0; i < 2; i++) {
            cl = CellRef(this, facecells[n][i]);
            cells[facecells[n][i]].faces[n] =
                new TreeFace(cl, cl.neighbour(n), n >> 1);
        }
    }
    #endif

    #ifdef TFTT_VERTICES
    // Add vertices

    // std::cout << "Adding Vertices.\n";

    // Corners
    TreeVertex* tv;
    for (int v = 0; v < 2*DIM; v++) {
        cells[v].vertices[v] = tv = new TreeVertex();
        tv->cells[~v & (2*DIM-1)] = CellRef(this, v);
        // std::cout << "\tcells[" << v << "].v[" << v << "] == " << v << "\n";
    }

    // Middle
    tv = new TreeVertex();
    for (int c = 0; c < 1<<DIM; c++) {
        tv->cells[c] = cl = CellRef(this, c);
        cells[c].vertices[~c & ((1<<DIM)-1)] = tv;

        // std::cout << "\tcells[" << cl << "].v[" << (~c & ((1<<DIM)-1)) << "] == " << c << "\n";
    }

    // Centers
    // Todo: Generalise

    // Left
    tv = new TreeVertex();
    tv->cells[1] = CellRef(this, 0);
    cells[0].vertices[2] = tv;
    tv->cells[3] = CellRef(this, 2);
    cells[2].vertices[0] = tv;

    // Right
    tv = new TreeVertex();
    tv->cells[2] = CellRef(this, 3);
    cells[3].vertices[1] = tv;
    tv->cells[0] = CellRef(this, 1);
    cells[1].vertices[3] = tv;

    // Top
    tv = new TreeVertex();
    tv->cells[0] = CellRef(this, 2);
    cells[2].vertices[3] = tv;
    tv->cells[1] = CellRef(this, 3);
    cells[3].vertices[2] = tv;

    // Bottom
    tv = new TreeVertex();
    tv->cells[2] = CellRef(this, 0);
    cells[0].vertices[1] = tv;
    tv->cells[3] = CellRef(this, 1);
    cells[1].vertices[0] = tv;

    #endif
}

TreeGroup::TreeGroup(int b)
    : parent()
{

    for (int i = 0; i < (1<<DIM); i++) {
        cells[i].children = nullptr;
        cells[i].rank = -1;
        cells[i].index = i;
    }

    id = ident_t::boundary(b);
    for (int i = 0; i < DIM; i++) {
        origin[i] = 0.0;
    }

    origin[b >> 1] = gtree.size[b >> 1];
    if ((b & 1) == 0)
        origin[b >> 1] *= -1;

    // Thread
    orientation = 0;
    next = CellRef();
    prev = CellRef();
}

TreeGroup::TreeGroup(CellRef p)
    : parent(p)
{
    #ifdef TFTT_DEBUG
    if (!p.isValid() || p.hasChildren()) {
        throw std::invalid_argument("Invalid cell ref for group parent.");
    }
    #endif


    for (int i = 0; i < (1<<DIM); i++) {
        cells[i].children = nullptr;
        cells[i].rank = p.rank();
        cells[i].index = i;
    }

    id = p.id().firstchild();

    for (int i = 0; i < DIM; i++) {
        origin[i] = p.origin(i);
    }

    // Update FTT
    if (p.group() == gtree.boundGroups) {
        for (int n = 0; n < 2*DIM; n++) {
            if (isBoundary() && n != (boundary() ^ 1)) continue;

            neighbours[n] = CellRef(-1);
        }
    }
    else {
        for (int n = 0; n < 2*DIM; n++) {
            if (isBoundary() && n != (boundary() ^ 1)) continue;

            neighbours[n] = p.neighbour(n);
            if (neighbours[n].hasChildren()) {
                neighbours[n].children()->neighbours[n ^ 1] = p;
            }
        }
    }

    // Thread
    if (!isBoundary()) {
        orientation = p.orientation();
        next = p.next();
        prev = p.prev();

        if (prev.isValid()) {
            if (prev.group()->next == p)
                prev.group()->next = CellRef(this, hilbChild(orientation, 0));
        }
        else {
            gtree.first = CellRef(this, hilbChild(orientation, 0));
        }

        if (next.isValid()) {
            if (next.group()->prev == p)
                next.group()->prev = CellRef(this, hilbChild(orientation, (1<<DIM)-1));
        }
        else {
            gtree.last = CellRef(this, hilbChild(orientation, (1<<DIM)-1));
        }

        if (gtree.firstActive == p) {
            gtree.firstActive = CellRef(this, hilbChild(orientation, 0));
        }
        if (gtree.lastActive == p) {
            gtree.lastActive = CellRef(this, hilbChild(orientation, (1<<DIM)-1));
        }
    }

    #ifdef TFTT_FACES
    // Update faces
    // Add internal faces

    // Todo: Generalise:
    // for (int n = 0; n < 4*(DIM-1); n++) {
    //     fc = new TreeFace();
    // }

    constexpr int pairs[4][4] = {
        {0, 1, 1, 0},
        {0, 2, 2, 3},
        {1, 2, 3, 3},
        {2, 1, 3, 0}
    };

    for (int i = 0; i < 4; i++) {
        cells[pairs[i][0]].faces[pairs[i][1]] =
            cells[pairs[i][2]].faces[pairs[i][3]] =
                new TreeFace(CellRef(this, pairs[i][0]), CellRef(this, pairs[i][2]), pairs[i][1] >> 1);
    }

    // Add external faces
    constexpr int facecells[4][2] = {
        {0, 2},
        {1, 3},
        {0, 1},
        {2, 3}
    };

    cell_t cl;
    for (int n = 0; n < 2*DIM; n++) {
        if (isBoundary() && n != (boundary() ^ 1)) continue;

        if (neighbours[n].level() <= p.level()) {
            // Add faces
            for (int i = 0; i < 2; i++) {
                cl = CellRef(this, facecells[n][i]);
                cells[facecells[n][i]].faces[n] =
                    new TreeFace(cl, cl.neighbour(n), n >> 1);
            }
        }
        else {
            // Re-use faces
            for (int i = 0; i < 2; i++) {
                cl = CellRef(this, facecells[n][i]);
                cells[facecells[n][i]].faces[n] = cl.neighbour(n).face(n ^ 1).face;
            }
        }
    }
    #endif

    if (!p.isBoundary()) {

        #ifdef TFTT_VERTICES
        // Add vertices
        // Todo: Generalise
        TreeVertex* tv;

        // Middle
        tv = new TreeVertex();
        for (int c = 0; c < 1<<DIM; c++) {
            tv->cells[c] = cl = CellRef(this, c);
            cells[c].vertices[~c & ((1<<DIM)-1)] = tv;
        }

        // Corners
        for (int v = 0; v < 1<<DIM; v++) {
            tv = p.vertex(v).vert;
            cells[v].vertices[v] = tv;
            tv->cells[~v & ((1<<DIM)-1)] = CellRef(this, v);
        }

        // Face Centers
        // constexpr int faceVertices[4][2] = {
        //     {1, 3},
        //     {0, 2},
        //     {2, 3},
        //     {0, 1}
        // };

        // CellRef nb;
        // for (int d = 0; d < 2*DIM; d++) {
        //     nb = p.neighbour(d);
        //     if (nb.hasChildren()) {
        //         // Already exists
        //         tv = nb.child(faceVertices[d][0]).vertex(faceVertices[d][1]).vert;
        //     }
        //     else {
        //         tv = new TreeVertex();
        //         for (int i = 0; i < 2*(DIM-1); i++) {
        //             tv->cells[CellRef::childIndexOnFace(d, 0)] = nb;
        //         }
        //     }

        //     for (int i = 0; i < 2*(DIM-1); i++) {
        //         tv->cells[CellRef::childIndexOnFace(~d & ((1<<DIM)-1), i)]
        //             = p.childOnFace(d, i);
        //     }
        // }



        // Left
        cell_t nb = p.neighbour(0);
        if (!nb.isBoundary() && nb.hasChildren()) {
            tv = nb.child(1).vertex(3).vert;
        }
        else {
            tv = new TreeVertex();
            tv->cells[0] = nb;
            tv->cells[2] = nb;
        }
        tv->cells[1] = CellRef(this, 0);
        cells[0].vertices[2] = tv;
        tv->cells[3] = CellRef(this, 2);
        cells[2].vertices[0] = tv;

        // Right
        nb = p.neighbour(1);
        if (!nb.isBoundary() && nb.hasChildren()) {
            tv = nb.child(0).vertex(2).vert;
        }
        else {
            tv = new TreeVertex();
            tv->cells[1] = nb;
            tv->cells[3] = nb;
        }
        tv->cells[2] = CellRef(this, 3);
        cells[3].vertices[1] = tv;
        tv->cells[0] = CellRef(this, 1);
        cells[1].vertices[3] = tv;

        // Bottom
        nb = p.neighbour(2);
        if (!nb.isBoundary() && nb.hasChildren()) {
            tv = nb.child(2).vertex(3).vert;
        }
        else {
            tv = new TreeVertex();
            tv->cells[0] = nb;
            tv->cells[1] = nb;
        }
        tv->cells[2] = CellRef(this, 0);
        cells[0].vertices[1] = tv;
        tv->cells[3] = CellRef(this, 1);
        cells[1].vertices[0] = tv;

        // Top
        nb = p.neighbour(3);
        if (!nb.isBoundary() && nb.hasChildren()) {
            tv = nb.child(0).vertex(1).vert;
        }
        else {
            tv = new TreeVertex();
            tv->cells[2] = nb;
            tv->cells[3] = nb;
        }
        tv->cells[0] = CellRef(this, 2);
        cells[2].vertices[3] = tv;
        tv->cells[1] = CellRef(this, 3);
        cells[3].vertices[2] = tv;

        #endif
    }
}

TreeGroup::~TreeGroup()
{
    #ifdef TFTT_FACES
    constexpr int facecells[4][2] = {
        {0, 2},
        {1, 3},
        {0, 1},
        {2, 3}
    };
    #endif

    if (!gtree.destroying && !isBoundary()) {
        if (prev.isValid()) {
            if (prev.group() != parent.group()) {
                if (prev.isLastInGroup()) {
                    prev.group()->next = parent;
                }
                if (parent.isFirstInGroup()) {
                    parent.group()->prev = prev;
                }
            }
        }
        else {
            gtree.first = parent;
        }

        if (next.isValid()) {
            if (next.group() != parent.group()) {
                if (next.isFirstInGroup()) {
                    next.group()->prev = parent;
                }
                if (parent.isLastInGroup()) {
                    parent.group()->next = next;
                }
            }
        }
        else {
            gtree.last = parent;
        }

        #ifdef TFTT_FACES
        // Faces
        for (int n = 0; n < 2*DIM; n++) {
            if (neighbours[n].level() < id.level()) {
                for (int i = 0; i < 2; i++) {
                    delete cells[facecells[n][i]].faces[n];
                }
            }
        }
        #endif
    }
    else {
        // delete faces nicely?
        // Should be ok to delete memory twice.
        // Todo: Make more efficient.

        // for (int n = 0; n < 2*DIM; n++) {
        //     for (int i = 0; i < 2*DIM; i++) {
        //         delete cells[n].faces[i];
        //     }
        // }

        // Delete vertices
        // delete cells[0].vertices[3];
    }

    for (unsigned int ch = 0; ch < cells.size(); ch++) {
        if (cells[ch].children)
            group::free(cells[ch].children);
    }
}


} // namespace tftt
