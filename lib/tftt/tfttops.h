
#ifndef TFTT_TFTTOPS_H
#define TFTT_TFTTOPS_H


#include "config.h"
#include "cellref.h"


namespace tftt {


cell_t find(ident_t id);
cell_t insert(ident_t id);
cell_t findmax(fnData dt, double* maxValRet);

void calcFaceCoefs(cell_t cl);

cell_t max(fnData dfn);
void relax(double omega, fnDataRef datafn, fnCell cellfn);
double resid(fnDataRef datafn, fnCell cellfn);


} // namespace tftt


#endif
