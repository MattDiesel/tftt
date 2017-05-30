
#ifdef TFTT_FACES

#include "../cellref.h"

#include "treeface.h"


namespace tftt {


TreeFace::TreeFace(cell_t ca, cell_t cb, int dim)
    : c1(ca), c2(cb), dimension(dim)
{
}


} // namespace tftt


#endif
