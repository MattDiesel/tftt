
#ifndef TFTT_TREECELL_H
#define TFTT_TREECELL_H


#include <array>

#include "tftt.h"


namespace tftt {

struct TreeGroup;
struct TreeFace;
struct TreeVertex;


struct TreeCell {
	data_t data;
    std::array<TreeFace*, DIM*2> faces;
    std::array<TreeVertex*, 1<<DIM> vertices;
    node_t rank;
	TreeGroup* children;

    int poisNgbC;
    CellRef poisNgb[12];
    bool poisNgbDir[12];
    double poisCoef[12];
    double poisAlpha[4];
};


} // namespace tftt


#endif
