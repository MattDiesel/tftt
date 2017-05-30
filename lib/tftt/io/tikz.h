
#ifndef TFTT_IO_TIKZ_H
#define TFTT_IO_TIKZ_H


#include <string>
#include <ostream>

#include "../config.h"
#include "../cellref.h"


namespace tftt {
namespace tikz {


void meshLayer(std::ostream& os, cell_t cl, int layer);
void meshLayer(std::ostream& ofs, int layer);
void meshLayer(std::string fname, int layer);

void mesh(std::string fname, int maxDepth);

void hilbert(std::string fname);

void morton(std::string fname);


} // namespace tikz
} // namespace tftt


#endif
