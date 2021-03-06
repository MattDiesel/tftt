
#include "../config.h"
#include "../cellref.h"
#include "../structure/tree.h"

#include "activecurve.h"


namespace tftt {


tagActiveCurve activecurve;


tagActiveCurve::curve_iterator tagActiveCurve::begin()
{
    return tagActiveCurve::curve_iterator(gtree.firstActive);
}


tagActiveCurve::curve_iterator tagActiveCurve::end()
{
    curve_iterator ret = tagActiveCurve::curve_iterator(gtree.lastActive);
    ret++;
    return ret;
}


tagActiveCurve::curve_iterator tagActiveCurve::rbegin()
{
    return tagActiveCurve::curve_iterator(gtree.lastActive);
}


tagActiveCurve::curve_iterator tagActiveCurve::rend()
{
    curve_iterator ret = tagActiveCurve::curve_iterator(gtree.firstActive);
    ret--;
    return ret;
}


} // namespace tftt
