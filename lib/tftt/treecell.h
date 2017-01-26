

#ifndef TFTT_TREECELL_H
#define TFTT_TREECELL_H


#include "tftt.h"


namespace tftt {

struct TreeGroup;


struct TreeCell {
	data_t data;
	TreeGroup* children;
};


} // namespace tftt


#endif
