
#ifndef TFTT_TREECELL_H
#define TFTT_TREECELL_H


#include <array>

#include "tftt.h"
#include "treeface.h"


namespace tftt {

struct TreeGroup;
struct TreeFace;


struct TreeCell {
	data_t data;
    std::array<TreeFace*, DIM*2> faces;
    node_t rank;
	TreeGroup* children;
};


} // namespace tftt


#endif
