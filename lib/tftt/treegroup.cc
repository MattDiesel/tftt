
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
            if (prev.group != p.group)
                prev.group->next = CellRef(this, hilbChild(orientation, 0));
        }
        else {
            gtree.first = CellRef(this, hilbChild(orientation, 0));
        }

        if (next.isValid()) {
            if (next.group != p.group)
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
}

TreeGroup::~TreeGroup() {
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
    }

    for (int ch = 0; ch < cells.size(); ch++) {
        if (cells[ch].children)
            delete cells[ch].children;
    }
}


} // namespace tftt
