
#ifndef TFTT_ITER_ACTIVECURVE_H
#define TFTT_ITER_ACTIVECURVE_H


#include "../config.h"
#include "curve.h"


namespace tftt {


struct tagActiveCurve : public tagCurve {
    curve_iterator begin();
    curve_iterator end();

    curve_iterator rbegin();
    curve_iterator rend();
};

extern tagActiveCurve activecurve;


} // namespace tftt


#endif
