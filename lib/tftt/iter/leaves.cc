
#include <iostream>

#include "../config.h"
#include "../cellref.h"
#include "../structure/tree.h"

#include "leaves.h"


namespace tftt {


tagLeaves leaves;

void tagLeaves::leaf_iterator::next()
{
    // if (!cr.isValid()) {
    //     return;
    // }

    if (cr.index()+1 >= 2*DIM) {
        if (cr.group() == gtree.root) {
            cr = CellRef();
            return;
        }
        cr = cr.parent();
        next();
        return;
    }
    else {
        cr.stepChild();
    }

    while (this->cr.hasChildren()) {
        this->cr = this->cr.child(0);
    }
}


tagLeaves::leaf_iterator::leaf_iterator(CellRef c) : cellref_iterator(c)
{
    if (!cr.isValid()) {
        return;
    }
    while (cr.hasChildren()) {
        cr = cr.child(0);
    }
}


tagLeaves::leaf_iterator tagLeaves::leaf_iterator::operator++()
{
    auto i = *this;
    next();
    return i;
}


tagLeaves::leaf_iterator tagLeaves::leaf_iterator::operator++(int junk)
{
    next();
    return *this;
}



tagLeaves::leaf_iterator tagLeaves::begin()
{
    return tagLeaves::leaf_iterator(CellRef(gtree.root, 0));
}


tagLeaves::leaf_iterator tagLeaves::end()
{
    return tagLeaves::leaf_iterator(CellRef());
}


} // namespace tftt

