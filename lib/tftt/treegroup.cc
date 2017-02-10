
#ifdef TFTT_DEBUG
	#include <stdexcept>
#endif

#include <iostream>

#include "tftt.h"
#include "cellref.h"
#include "tree.h"

#include "treegroup.h"


namespace tftt {

TreeGroup::TreeGroup()
		: parent(), boundary(-1) {
	for (int i = 0; i < (1<<DIM); i++) {
		cells[i].children = nullptr;
	}

	id = 0;
	for (int i = 0; i < DIM; i++) {
		origin[i] = 0.0;
	}

	for (int i = 0; i < DIM*2; i++) {
		neighbours[i] = CellRef(gtree.boundGroups, i);
	}

	flaggedForCoarsening = false;
}

TreeGroup::TreeGroup(int b)
		: parent(), boundary(b) {

	for (int i = 0; i < (1<<DIM); i++) {
		cells[i].children = nullptr;
	}

	id = 0;
	for (int i = 0; i < DIM; i++) {
		origin[i] = 0.0;
	}

	origin[b >> 1] = gtree.size[b >> 1];
	if (b & 1 == 0)
		origin[b >> 1] *= -1;

	flaggedForCoarsening = false;
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
}

TreeGroup::~TreeGroup() {
	for (int ch = 0; ch < cells.size(); ch++) {
		if (cells[ch].children)
			delete cells[ch].children;
	}
}


} // namespace tftt
