
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


enum ADAPTFLAGS {
    AF_NoAction = 0,
    AF_Refine,
    AF_Coarsen,
    AF_HoldRefined,
    AF_HoldCoarsened
};

struct TreeCell {
    data_t data;
    node_t rank;
    TreeGroup* children;

    #ifdef TFTT_FACES
    std::array<TreeFace*, DIM*2> faces;
    #endif

    #ifdef TFTT_VERTICES
    std::array<TreeVertex*, 1<<DIM> vertices;
    #endif

    int poisNgbC;
    CellRef poisNgb[12];
    double poisCoef[12];
    double poisAlpha[4];
    double cenCoef;

    ADAPTFLAGS adaptFlags;
    uint8_t adaptVector[2*DIM];
    uint8_t adaptHoldVector[2*DIM];
};


} // namespace tftt


#endif
