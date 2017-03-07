
#ifndef TFTT_TREEFACE_H
#define TFTT_TREEFACE_H


#include "tftt.h"


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


