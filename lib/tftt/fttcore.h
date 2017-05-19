
#ifndef TFTT_FTTCORE_H
#define TFTT_FTTCORE_H


#include "config.h"
#include "cellref.h"


namespace tftt {


//! Constructs the top level tree
void init(double w, double h);
void reset();

void refine(CellRef cl);
void coarsen(CellRef cl);

CellRef atPos(double pos[DIM]);
CellRef atVertex(int v);


} // namespace tftt


#endif
