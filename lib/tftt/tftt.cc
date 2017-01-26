
#include <iostream> // TODO: Remove
#include <stdexcept>

#include "tftt.h"
#include "tree.h"
#include "cellref.h"


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

    // Neighbours
    gtree.root->neighbours[0] = &CopyNeighbour;
    gtree.root->neighbours[1] = &CopyNeighbour;
    for (int i = 2; i < 2*DIM; i++) {
    	gtree.root->neighbours[i] = &ReflectNeighbour;
    }
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


} // namespace tftt
