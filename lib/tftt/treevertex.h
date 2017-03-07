
#ifndef TFTT_TREEVERTEX_H
#define TFTT_TREEVERTEX_H


#include <array>

#include "tftt.h"


namespace tftt {


struct TreeVertex {
    vertexdata_t data;
    std::array<CellRef,1<<DIM> cells;
};


} // namespace tftt


#endif