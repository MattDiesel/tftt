
#include "../config.h"
#include "../cellref.h"
#include "../structure/tree.h"

#include "leaforthos.h"


namespace tftt {


tagLeafOrthos leaforthos;


void nextLeaf(cell_t& cr)
{
    if (cr.index()+1 >= 2*DIM) {
        if (cr.group() == gtree.root) {
            cr = CellRef();
            return;
        }
        cr = cr.parent();
        nextLeaf(cr);
        return;
    }
    else {
        cr.stepChild();
    }

    while (cr.hasChildren()) {
        cr = cr.child(0);
    }
}


void tagLeafOrthos::ortho_iterator::next()
{
    cell_t par = cr.parent();
    do {
        nextLeaf(cr);
    }
    while (cr.isValid() && cr.parent() == par);

    if (cr.isValid())
        ortho = cr.parent();
    else
        ortho = CellRef();
}


tagLeafOrthos::ortho_iterator::ortho_iterator(CellRef c) : cr(c)
{
    if (!cr.isValid()) {
        return;
    }
    while (cr.hasChildren()) {
        cr = cr.child(0);
    }
    ortho = cr.parent();
}


tagLeafOrthos::ortho_iterator tagLeafOrthos::ortho_iterator::operator++()
{
    auto i = *this;
    next();
    return i;
}


tagLeafOrthos::ortho_iterator tagLeafOrthos::ortho_iterator::operator++(int junk)
{
    next();
    return *this;
}


CellRef& tagLeafOrthos::ortho_iterator::operator*()
{
    return ortho;
}


CellRef* tagLeafOrthos::ortho_iterator::operator->()
{
    return &ortho;
}


bool tagLeafOrthos::ortho_iterator::operator==(const tagLeafOrthos::ortho_iterator& rhs)
{
    return ortho == rhs.ortho;
}


bool tagLeafOrthos::ortho_iterator::operator!=(const tagLeafOrthos::ortho_iterator& rhs)
{
    return !(ortho == rhs.ortho);
}


tagLeafOrthos::ortho_iterator tagLeafOrthos::begin()
{
    return tagLeafOrthos::ortho_iterator(CellRef(gtree.root, 0));
}


tagLeafOrthos::ortho_iterator tagLeafOrthos::end()
{
    return tagLeafOrthos::ortho_iterator(CellRef());
}


} // namespace tftt
