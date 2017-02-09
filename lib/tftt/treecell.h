
#ifndef TFTT_TREECELL_H
#define TFTT_TREECELL_H


#include <array>

#include "tftt.h"


namespace tftt {

struct TreeGroup;


struct TreeCell {
	data_t data;
    std::array<facedata_t, DIM*2> facedata;
	TreeGroup* children;
};


} // namespace tftt


#endif
