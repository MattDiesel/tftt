
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

void meshLayer_if(std::ostream& os, cell_t cl, int layer, fnCheckCell pred);
void meshLayer_if(std::ostream& ofs, int layer, fnCheckCell pred);
void meshLayer_if(std::string fname, int layer, fnCheckCell pred);

void mesh_if(std::string fname, int maxDepth, fnCheckCell pred);


} // namespace tikz
} // namespace tftt


#endif
