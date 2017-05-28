
#ifndef TFTT_TREECELL_H
#define TFTT_TREECELL_H


#include <array>

#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct TreeGroup;

#ifdef TFTT_FACES
    struct TreeFace;
#endif
#ifdef TFTT_VERTICES
    struct TreeVertex;
#endif


enum ADAPTFLAGS : uint8_t {
    AF_NoAction = 0,
    AF_Refine,
    AF_Coarsen,
    AF_HoldRefined,
    AF_HoldCoarsened
};

struct TreeCell {
    int8_t index;
    node_t rank;

    int8_t poisNgbC;

    ADAPTFLAGS adaptFlags;
    uint8_t adaptVector[2*DIM];
    uint8_t adaptHoldVector[2*DIM];

    data_t data;
    TreeGroup* children;

    #ifdef TFTT_FACES
    std::array<TreeFace*, DIM*2> faces;
    #endif

    #ifdef TFTT_VERTICES
    std::array<TreeVertex*, 1<<DIM> vertices;
    #endif

    CellRef poisNgb[8];
    double poisCoef[8];
    double cenCoef;

};


} // namespace tftt


#endif
