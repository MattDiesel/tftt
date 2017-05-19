
#include "../config.h"
#include "../cellref.h"
#include "../structure/tree.h"

#include "neighbours.h"


namespace tftt {


tagNeighbours::tagNeighbours(CellRef c)
{
    cr = c;
}


void tagNeighbours::neighb_iterator::next()
{
    nb++;
    cr = cl.neighbour(nb);
}


tagNeighbours::neighb_iterator::neighb_iterator(CellRef c, int n)
{
    cl = c;
    nb = n;
    cr = c.neighbour(n);
}


tagNeighbours::neighb_iterator tagNeighbours::neighb_iterator::operator++()
{
    auto i = *this;
    next();
    return i;
}


tagNeighbours::neighb_iterator tagNeighbours::neighb_iterator::operator++(int junk)
{
    next();
    return *this;
}


CellRef& tagNeighbours::neighb_iterator::operator*()
{
    return cr;
}


CellRef* tagNeighbours::neighb_iterator::operator->()
{
    return &cr;
}


bool tagNeighbours::neighb_iterator::operator==(const tagNeighbours::neighb_iterator& rhs)
{
    return (cl == rhs.cl) && (nb == rhs.nb);
}


bool tagNeighbours::neighb_iterator::operator!=(const tagNeighbours::neighb_iterator& rhs)
{
    return (cl != rhs.cl) || (nb != rhs.nb);
}


tagNeighbours::neighb_iterator tagNeighbours::begin()
{
    return tagNeighbours::neighb_iterator(cr, 0);
}


tagNeighbours::neighb_iterator tagNeighbours::end()
{
    return tagNeighbours::neighb_iterator(cr, 2*DIM);
}


} // namespace tftt
