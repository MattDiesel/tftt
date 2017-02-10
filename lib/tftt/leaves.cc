
#include <iostream>

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
        return;
    }
    else {
        cr.index++;
    }

    while (this->cr.hasChildren()) {
        this->cr = this->cr.child(0);
    }
}

tagLeaves::leaf_iterator::leaf_iterator(CellRef c) : cr(c) {
    if (!cr.isValid()) {
        return;
    }
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





tagLeafOrthos leaforthos;



void nextLeaf(cell_t& cr) {
    if (cr.index+1 >= 2*DIM) {
        if (cr.group == gtree.root) {
            cr = CellRef();
            return;
        }
        cr = cr.parent();
        nextLeaf(cr);
        return;
    }
    else {
        cr.index++;
    }

    while (cr.hasChildren()) {
        cr = cr.child(0);
    }
}

void tagLeafOrthos::ortho_iterator::next() {
    cell_t par = cr.parent();
    do {
        nextLeaf(cr);
    } while (cr.isValid() && cr.parent() == par);

    if (cr.isValid())
        ortho = cr.parent();
    else
        ortho = CellRef();
}

tagLeafOrthos::ortho_iterator::ortho_iterator(CellRef c) : cr(c) {
    if (!cr.isValid()) {
        return;
    }
    while (cr.hasChildren()) {
        cr = cr.child(0);
    }
    ortho = cr.parent();
}

tagLeafOrthos::ortho_iterator tagLeafOrthos::ortho_iterator::operator++() {
    auto i = *this;
    next();
    return i;
}

tagLeafOrthos::ortho_iterator tagLeafOrthos::ortho_iterator::operator++(int junk) {
    next();
    return *this;
}

CellRef& tagLeafOrthos::ortho_iterator::operator*() {
    return ortho;
}

CellRef* tagLeafOrthos::ortho_iterator::operator->() {
    return &ortho;
}

bool tagLeafOrthos::ortho_iterator::operator==(const tagLeafOrthos::ortho_iterator& rhs) {
    return ortho == rhs.ortho;
}

bool tagLeafOrthos::ortho_iterator::operator!=(const tagLeafOrthos::ortho_iterator& rhs) {
    return !(ortho == rhs.ortho);
}

tagLeafOrthos::ortho_iterator tagLeafOrthos::begin() {
    return tagLeafOrthos::ortho_iterator(CellRef(gtree.root, 0));
}

tagLeafOrthos::ortho_iterator tagLeafOrthos::end() {
    return tagLeafOrthos::ortho_iterator(CellRef());
}





tagBoundaryLeaves boundaryCells(int b) {
    tagBoundaryLeaves ret;
    ret.b = b;
    return ret;
}


constexpr int _start[4] = {1,0,2,0};
constexpr int _inc[4]   = {2,2,1,1};
constexpr int _end[4]   = {3,2,3,1};

void tagBoundaryLeaves::bleaf_iterator::next() {
    if (!cr.isValid()) {
        return;
    }

    if (cr.index+_inc[b] > _end[b]) {
        if (cr.parent().group == gtree.boundGroups) {
            cr = CellRef();
            return;
        }
        cr = cr.parent();
        next();
        return;
    }
    else {
        cr.index += _inc[b];
    }

    while (this->cr.hasChildren()) {
        this->cr = this->cr.child(_start[b]);
    }
}

tagBoundaryLeaves::bleaf_iterator::bleaf_iterator(CellRef c, int bnd) : cr(c), b(bnd) {
    if (!cr.isValid()) {
        return;
    }
    while (cr.hasChildren()) {
        cr = cr.child(_start[b]);
    }
}

tagBoundaryLeaves::bleaf_iterator tagBoundaryLeaves::bleaf_iterator::operator++() {
    auto i = *this;
    next();
    return i;
}

tagBoundaryLeaves::bleaf_iterator tagBoundaryLeaves::bleaf_iterator::operator++(int junk) {
    next();
    return *this;
}

CellRef& tagBoundaryLeaves::bleaf_iterator::operator*() {
    return cr;
}

CellRef* tagBoundaryLeaves::bleaf_iterator::operator->() {
    return &cr;
}

bool tagBoundaryLeaves::bleaf_iterator::operator==(const tagBoundaryLeaves::bleaf_iterator& rhs) {
    return cr == rhs.cr;
}

bool tagBoundaryLeaves::bleaf_iterator::operator!=(const tagBoundaryLeaves::bleaf_iterator& rhs) {
    return !(cr == rhs.cr);
}

tagBoundaryLeaves::bleaf_iterator tagBoundaryLeaves::begin() {
    return tagBoundaryLeaves::bleaf_iterator(CellRef(gtree.boundGroups, b), b);
}

tagBoundaryLeaves::bleaf_iterator tagBoundaryLeaves::end() {
    return tagBoundaryLeaves::bleaf_iterator(CellRef(), b);
}




} // namespace tftt
