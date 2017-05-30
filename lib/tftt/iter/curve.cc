
#include "../config.h"
#include "../cellref.h"
#include "../structure/tree.h"

#include "curve.h"


namespace tftt {


tagCurve curve;


void tagCurve::curve_iterator::next()
{
    cr = cr.next();
}


void tagCurve::curve_iterator::prev()
{
    cr = cr.prev();
}


tagCurve::curve_iterator::curve_iterator(cell_t c) : cellref_iterator(c)
{
}


tagCurve::curve_iterator tagCurve::curve_iterator::operator++()
{
    auto i = *this;
    next();
    return i;
}


tagCurve::curve_iterator tagCurve::curve_iterator::operator++(int junk)
{
    next();
    return *this;
}


tagCurve::curve_iterator tagCurve::curve_iterator::operator--()
{
    auto i = *this;
    prev();
    return i;
}


tagCurve::curve_iterator tagCurve::curve_iterator::operator--(int junk)
{
    prev();
    return *this;
}


tagCurve::curve_iterator tagCurve::begin()
{
    return tagCurve::curve_iterator(gtree.first);
}


tagCurve::curve_iterator tagCurve::end()
{
    return tagCurve::curve_iterator(cell_t());
}


tagCurve::curve_iterator tagCurve::rbegin()
{
    return tagCurve::curve_iterator(gtree.last);
}


tagCurve::curve_iterator tagCurve::rend()
{
    return tagCurve::curve_iterator(cell_t());
}


} // namespace tftt
