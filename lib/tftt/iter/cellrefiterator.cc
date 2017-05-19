
#include "../config.h"
#include "../cellref.h"

#include "cellrefiterator.h"


namespace tftt {


cellref_iterator::cellref_iterator(CellRef c)
    : cr(c)
{
}


CellRef& cellref_iterator::operator*()
{
    return cr;
}


CellRef* cellref_iterator::operator->()
{
    return &cr;
}


bool cellref_iterator::operator==(const cellref_iterator& rhs)
{
    return cr == rhs.cr;
}


bool cellref_iterator::operator!=(const cellref_iterator& rhs)
{
    return !(cr == rhs.cr);
}


} // namespace tftt
