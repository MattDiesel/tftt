

#ifndef RT_COMPUTEPOISCEOF_H
#define RT_COMPUTEPOISCEOF_H


#include "tftt/tftt.h"


double faceVoF(tftt::cell_t const& c1, int dir);
double faceRho(tftt::cell_t const& c1, int dir);

void computePoisCoef();


#endif
