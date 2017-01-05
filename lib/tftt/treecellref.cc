
#ifdef TFTT_DEBUG
	#include <stdexcept>
#endif

#include "tftt/config.h"
#include "treecellref.h"


namespace tftt {


TreeCellRef::TreeCellRef(TreeGroup* gr, int ind)
		: group(gr), index(ind){
	#ifdef TFTT_DEBUG
		if (!gr) throw std::invalid_argument("Null TreeGroup given.");
		if (ind < 0 || ind > (1 << DIM)) throw std::invalid_argument("Child index is invalid for dimension of problem.");
	#endif
}


TreeCellRef::TreeCellRef(bool copy)
		: group(nullptr), index(0) {
	if (copy)
		index = TCR_COPYBOUNDARY;
}


bool TreeCellRef::isBoundary() const {
	return !group;
}


bool TreeCellRef::isCopyBoundary() const {
	return isBoundary() && (index & TCR_COPYBOUNDARY);
}


bool TreeCellRef::isReflectBoundary() const {
	return isBoundary() && !(index & TCR_COPYBOUNDARY);
}


} // namespace tftt
