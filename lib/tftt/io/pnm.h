
#ifndef TFTT_IO_PNM_H
#define TFTT_IO_PNM_H


#include <string>
#include <ostream>

#include "../config.h"
#include "../cellref.h"


namespace tftt {
namespace pnm {


//! Draws the normalised data to an image file
void pgm(std::string fname, int imgW, int imgH, fnCell data);
//! \copydoc pgm
void pgm(std::ostream& os, int imgW, int imgH, fnCell data);


} // namespace pnm
} // namespace tftt


#endif
