
#ifndef TFTT_TREEFACE_H
#define TFTT_TREEFACE_H

#ifdef TFTT_FACES

#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct TreeFace {
    facedata_t data;

    cell_t c1;
    cell_t c2;
    int dimension;

    TreeFace(cell_t ca, cell_t cb, int dim);
};


} // namespace tftt


#endif

#endif
