
#ifndef TFTT_TREEFACE_H
#define TFTT_TREEFACE_H

#ifdef TFTT_FACES

#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct TreeFace {
    facedata_t data;

    CellRef c1;
    CellRef c2;
    int dimension;

    TreeFace(CellRef ca, CellRef cb, int dim);
};


} // namespace tftt


#endif

#endif
