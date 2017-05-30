
#ifndef TFTT_FTTCORE_H
#define TFTT_FTTCORE_H


#include "config.h"
#include "cellref.h"


namespace tftt {


//! Constructs the top level tree
void init(double w, double h);
void reset();

void refine(cell_t cl);
void coarsen(cell_t cl);

cell_t atPos(double pos[DIM]);
cell_t atVertex(int v);

cell_t find(ident_t id);
cell_t insert(ident_t id);
cell_t findmax(fnData dt, double* maxValRet);
cell_t max(fnData dfn);


} // namespace tftt


#endif
