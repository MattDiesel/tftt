
#include "tftt.h"

#include "treeface.h"


namespace tftt {


TreeFace::TreeFace(CellRef ca, CellRef cb, int dim)
        : c1(ca), c2(cb), dimension(dim) {
}


} // namespace tftt
