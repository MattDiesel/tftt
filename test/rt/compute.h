
#ifndef RT_COMPUTE_H
#define RT_COMPUTE_H

#include "computePoisCoef.h"

void computeRho();
void computeAdv();
void computeVisc();
void computeIntermVelo();
void computeFlux();
void computeDive();
double computeCompatability();

void correctVelo();
void correctFlux();


#endif
