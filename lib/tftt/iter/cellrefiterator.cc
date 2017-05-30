
#include "../config.h"
#include "../cellref.h"

#include "cellrefiterator.h"


namespace tftt {


cellref_iterator::cellref_iterator(cell_t c)
    : cr(c)
{
}


cell_t& cellref_iterator::operator*()
{
    return cr;
}


cell_t* cellref_iterator::operator->()
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
