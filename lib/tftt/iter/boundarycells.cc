
#include "../config.h"
#include "../cellref.h"
#include "../structure/tree.h"

#include "boundarycells.h"


namespace tftt {


tagBoundaryLeaves boundaryCells(int b)
{
    tagBoundaryLeaves ret;
    ret.b = b;
    return ret;
}


constexpr int _start[4] = {1,0,2,0};
constexpr int _inc[4]   = {2,2,1,1};
constexpr int _end[4]   = {3,2,3,1};

void tagBoundaryLeaves::bleaf_iterator::next()
{
    if (!cr.isValid()) {
        return;
    }

    if (cr.index()+_inc[b] > _end[b]) {
        if (cr.parent().group() == gtree.boundGroups) {
            cr = cell_t();
            return;
        }
        cr = cr.parent();
        next();
        return;
    }
    else {
        cr = cell_t(cr.group(), cr.index()+_inc[b]);
    }

    while (this->cr.hasChildren()) {
        this->cr = this->cr.child(_start[b]);
    }
}


tagBoundaryLeaves::bleaf_iterator::bleaf_iterator(cell_t c, int bnd)
    : cellref_iterator(c), b(bnd)
{
    if (!cr.isValid()) {
        return;
    }
    while (cr.hasChildren()) {
        cr = cr.child(_start[b]);
    }
}


tagBoundaryLeaves::bleaf_iterator tagBoundaryLeaves::bleaf_iterator::operator++()
{
    auto i = *this;
    next();
    return i;
}


tagBoundaryLeaves::bleaf_iterator tagBoundaryLeaves::bleaf_iterator::operator++(int junk)
{
    next();
    return *this;
}


tagBoundaryLeaves::bleaf_iterator tagBoundaryLeaves::begin()
{
    return tagBoundaryLeaves::bleaf_iterator(cell_t(gtree.boundGroups, b), b);
}


tagBoundaryLeaves::bleaf_iterator tagBoundaryLeaves::end()
{
    return tagBoundaryLeaves::bleaf_iterator(cell_t(), b);
}


} // namespace tftt
