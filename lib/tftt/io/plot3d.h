
#ifndef TFTT_IO_PLOT3D_H
#define TFTT_IO_PLOT3D_H


#include <string>
#include <ostream>

#include "../config.h"
#include "../cellref.h"


namespace tftt {
namespace plot3d {


//! Plots the data
void scatter(std::string fname, fnCell dt);
//! \copydoc scatter
void scatter(std::ostream& os, fnCell dt);


//! Plots the data
void sample(std::string fname, int resX, int resY, fnCell dt);
//! \copydoc sample
void sample(std::ostream& os, int resX, int resY, fnCell dt);


} // namespace plot3d
} // namespace tftt


#endif
