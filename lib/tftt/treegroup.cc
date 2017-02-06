
#ifdef TFTT_DEBUG
	#include <stdexcept>
#endif

#include "tftt.h"
#include "cellref.h"

#include "treegroup.h"


namespace tftt {

TreeGroup CopyNeighbour;
TreeGroup ReflectNeighbour;


TreeGroup::TreeGroup()
		: parent() {
	for (int i = 0; i < (1<<DIM); i++) {
		cells[i].children = nullptr;
	}

	id = 0;
	for (int i = 0; i < DIM; i++) {
		origin[i] = 0.0;
	}

	neighbours[0] = CellRef(true);
	neighbours[1] = CellRef(true);
	for (int i = 2; i < DIM*2; i++) {
		neighbours[i] = CellRef(false);
	}

	flaggedForCoarsening = false;
}


TreeGroup::TreeGroup(CellRef p)
		: parent(p) {

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
	}
}

TreeGroup::~TreeGroup() {
}


bool TreeGroup::isCopyBoundary() const {
	return this == &CopyNeighbour;
}


bool TreeGroup::isReflectBoundary() const {
	return this == &ReflectNeighbour;
}


bool TreeGroup::isBoundary() const {
	return isCopyBoundary() || isReflectBoundary();
}





} // namespace tftt
