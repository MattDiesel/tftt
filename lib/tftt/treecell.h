

#ifndef TFTT_TREECELL_H
#define TFTT_TREECELL_H


#include "tftt/config.h"
#include "treegroup.h"


namespace tftt {


struct TreeCell {
	data_t data;
	TreeGroup* children;
};


} // namespace tftt


#endif
