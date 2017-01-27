
#include "tftt.h"
#include "tree.h"

#include "leaves.h"

namespace tftt {


tagLeaves leaves;

void tagLeaves::leaf_iterator::next() {
    // if (!cr.isValid()) {
    //     return;
    // }

    if (cr.index+1 >= 2*DIM) {
        if (cr.group == gtree.root) {
            cr = CellRef();
            return;
        }
        cr = cr.parent();
        next();
    }
    else {
        cr.index++;
    }

    while (cr.hasChildren()) {
        cr = cr.child(0);
    }
}

tagLeaves::leaf_iterator::leaf_iterator(CellRef c) : cr(c) {
    while (cr.hasChildren()) {
        cr = cr.child(0);
    }
}

tagLeaves::leaf_iterator tagLeaves::leaf_iterator::operator++() {
    auto i = *this;
    next();
    return i;
}

tagLeaves::leaf_iterator tagLeaves::leaf_iterator::operator++(int junk) {
    next();
    return *this;
}

CellRef& tagLeaves::leaf_iterator::operator*() {
    return cr;
}

CellRef* tagLeaves::leaf_iterator::operator->() {
    return &cr;
}

bool tagLeaves::leaf_iterator::operator==(const tagLeaves::leaf_iterator& rhs) {
    return cr == rhs.cr;
}

bool tagLeaves::leaf_iterator::operator!=(const tagLeaves::leaf_iterator& rhs) {
    return !(cr == rhs.cr);
}

tagLeaves::leaf_iterator tagLeaves::begin() {
    return tagLeaves::leaf_iterator(CellRef(gtree.root, 0));
}

tagLeaves::leaf_iterator tagLeaves::end() {
    return tagLeaves::leaf_iterator(CellRef());
}


} // namespace tftt
