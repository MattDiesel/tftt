
#ifdef TFTT_DEBUG
    #include <stdexcept>
#endif

#include <iostream>

#include "tftt.h"
#include "cellref.h"
#include "tree.h"
#include "hilbert.h"

#include "treegroup.h"


namespace tftt {

TreeGroup::TreeGroup()
        : parent(), boundary(-1) {
    for (int i = 0; i < (1<<DIM); i++) {
        cells[i].children = nullptr;
        cells[i].rank = -1;
    }

    id = 0;
    for (int i = 0; i < DIM; i++) {
        origin[i] = 0.0;
    }

    for (int i = 0; i < DIM*2; i++) {
        neighbours[i] = CellRef(gtree.boundGroups, i);
    }

    flaggedForCoarsening = false;

    // Thread
    orientation = 0;
    next = CellRef();
    prev = CellRef();

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
}

TreeGroup::TreeGroup(int b)
        : parent(), boundary(b) {

    for (int i = 0; i < (1<<DIM); i++) {
        cells[i].children = nullptr;
        cells[i].rank = -1;
    }

    id = 0;
    for (int i = 0; i < DIM; i++) {
        origin[i] = 0.0;
    }

    origin[b >> 1] = gtree.size[b >> 1];
    if (b & 1 == 0)
        origin[b >> 1] *= -1;

    flaggedForCoarsening = false;

    // Thread
    orientation = 0;
    next = CellRef();
    prev = CellRef();
}

TreeGroup::TreeGroup(CellRef p)
        : parent(p) {

    boundary = p.group->boundary;

    #ifdef TFTT_DEBUG
        if (!p.isValid() || p.hasChildren()) {
            throw std::argument_exception("Invalid cell ref for group parent.");
        }
    #endif


    for (int i = 0; i < (1<<DIM); i++) {
        cells[i].children = nullptr;
        cells[i].rank = p.rank();
    }

    id = p.id().firstchild();

    for (int i = 0; i < DIM; i++) {
        origin[i] = p.origin(i);
    }

    // Update FTT
    for (int n = 0; n < 2*DIM; n++) {
        neighbours[n] = p.neighbour(n);
        if (neighbours[n].hasChildren()) {
            neighbours[n].children()->neighbours[n ^ 1] = p;
        }
    }

    flaggedForCoarsening = false;

    // Thread
    if (boundary == -1) {
        orientation = p.orientation();
        next = p.next();
        prev = p.prev();

        if (prev.isValid()) {
            if (prev.group->next == p)
                prev.group->next = CellRef(this, hilbChild(orientation, 0));
        }
        else {
            gtree.first = CellRef(this, hilbChild(orientation, 0));
        }

        if (next.isValid()) {
            if (next.group->prev == p)
                next.group->prev = CellRef(this, hilbChild(orientation, (1<<DIM)-1));
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
}

TreeGroup::~TreeGroup() {
    constexpr int facecells[4][2] = {
        {0, 2},
        {1, 3},
        {0, 1},
        {2, 3}
    };

    if (!gtree.destroying && boundary == -1) {
        if (prev.isValid()) {
            if (prev.group != parent.group) {
                if (prev.isLastInGroup()) {
                    prev.group->next = parent;
                }
                if (parent.isFirstInGroup()) {
                    parent.group->prev = prev;
                }
            }
        }
        else {
            gtree.first = parent;
        }

        if (next.isValid()) {
            if (next.group != parent.group) {
                if (next.isFirstInGroup()) {
                    next.group->prev = parent;
                }
                if (parent.isLastInGroup()) {
                    parent.group->next = next;
                }
            }
        }
        else {
            gtree.last = parent;
        }

        // Faces 
        for (int n = 0; n < 2*DIM; n++) {
            if (neighbours[n].level() < id.level()) {
                for (int i = 0; i < 2; i++) {
                    delete cells[facecells[n][i]].faces[n];
                }
            }
        }
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
    }

    for (int ch = 0; ch < cells.size(); ch++) {
        if (cells[ch].children)
            delete cells[ch].children;
    }
}


} // namespace tftt
