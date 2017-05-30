
#ifndef TFTT_TREEVERTEX_H
#define TFTT_TREEVERTEX_H

#ifdef TFTT_VERTICES


#include <array>

#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct TreeVertex {
    vertexdata_t data;
    std::array<cell_t,1<<DIM> cells;
};


} // namespace tftt


#endif
#endif
