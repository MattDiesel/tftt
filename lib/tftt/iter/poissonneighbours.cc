
#include "../config.h"
#include "../cellref.h"
#include "../structure/tree.h"

#include "poissonneighbours.h"


namespace tftt {


tagPoissonNeighbours::tagPoissonNeighbours(CellRef c)
{
    cr = c;
}


void tagPoissonNeighbours::poisneighb_iterator::next()
{
    nb++;

    TreeCell& tc = cl.group()->cells[cl.index()];
    cr = tc.poisNgb[nb];
}


tagPoissonNeighbours::poisneighb_iterator::poisneighb_iterator(CellRef c, int n)
{
    cl = c;
    nb = n;

    TreeCell& tc = cl.group()->cells[cl.index()];
    cr = tc.poisNgb[nb];
}


tagPoissonNeighbours::poisneighb_iterator tagPoissonNeighbours::poisneighb_iterator::operator++()
{
    auto i = *this;
    next();
    return i;
}


tagPoissonNeighbours::poisneighb_iterator tagPoissonNeighbours::poisneighb_iterator::operator++(int junk)
{
    next();
    return *this;
}


CellRef& tagPoissonNeighbours::poisneighb_iterator::operator*()
{
    return cr;
}


CellRef* tagPoissonNeighbours::poisneighb_iterator::operator->()
{
    return &cr;
}


bool tagPoissonNeighbours::poisneighb_iterator::operator==(const tagPoissonNeighbours::poisneighb_iterator& rhs)
{
    return (cl == rhs.cl) && (nb == rhs.nb);
}


bool tagPoissonNeighbours::poisneighb_iterator::operator!=(const tagPoissonNeighbours::poisneighb_iterator& rhs)
{
    return (cl != rhs.cl) || (nb != rhs.nb);
}


tagPoissonNeighbours::poisneighb_iterator tagPoissonNeighbours::begin()
{
    return tagPoissonNeighbours::poisneighb_iterator(cr, 0);
}


tagPoissonNeighbours::poisneighb_iterator tagPoissonNeighbours::end()
{
    TreeCell& tc = cr.group()->cells[cr.index()];
    return tagPoissonNeighbours::poisneighb_iterator(cr, tc.poisNgbC);
}


} // namespace tftt
