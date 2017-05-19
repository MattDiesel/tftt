
#ifndef TFTT_TFTTOPS_H
#define TFTT_TFTTOPS_H


#include "config.h"
#include "cellref.h"


namespace tftt {

void calcFaceCoefs(cell_t cl);

void relax(double omega, fnDataRef datafn, fnCell cellfn);
double resid(fnDataRef datafn, fnCell cellfn);


} // namespace tftt


#endif
