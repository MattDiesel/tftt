
#ifdef TFTT_DEBUG
    #include <stdexcept>
#endif

#include <ostream>

#include "tree.h"
#include "treecell.h"

#include "cellref.h"


namespace tftt {


CellRef::CellRef()
        : group(nullptr), index(0) {
}


CellRef::CellRef(TreeGroup* gr, int ind)
        : group(gr), index(ind){
    #ifdef TFTT_DEBUG
        if (!gr) throw std::invalid_argument("Null TreeGroup given.");
        if (ind < 0 || ind > (1 << DIM)) throw std::invalid_argument("Child index is invalid for dimension of problem.");
    #endif
}


CellRef::CellRef(bool copy)
        : group(nullptr), index(0) {
    if (copy)
        index = TCR_COPYBOUNDARY;
    else
        index = TCR_REFLECTBOUNDARY;
}


bool CellRef::isValid() const {
    return group || index;
}


bool CellRef::isBoundary() const {
    return !group && index;
}


bool CellRef::isCopyBoundary() const {
    return isBoundary() && (index & TCR_COPYBOUNDARY);
}


bool CellRef::isReflectBoundary() const {
    return isBoundary() && (index & TCR_REFLECTBOUNDARY);
}

bool CellRef::flaggedForCoarsening() const {
    return group->flaggedForCoarsening;
}


CellRef CellRef::parent() const {
    return group->parent;
}

CellRef CellRef::child(int n) const {
    return CellRef(children(), n);
}

bool CellRef::hasChildren() const {
    return group && children();
}

TreeGroup* CellRef::children() const {
    return group->cells[index].children;
}


bool nbInParent(int ch, int nb) {
    // Explanation:
    // Reduce to 1d problem, then check bounds.
    
    return ((ch >> (nb >> 1)) ^ nb) & 1;
}


CellRef CellRef::neighbour(int n) const {
    if (nbInParent(index, n))
        return CellRef(group, index ^ (1 << (n >> 1)));

    CellRef cr = group->neighbours[n];
    if (cr.isBoundary() || !cr.hasChildren())
        return cr;

    return cr.child(index ^ (1 << (n >> 1)));
}


ident_t CellRef::id() const {
    return group->id.id | index;
}

int CellRef::level() const {
    return group->id.level();
}


data_t& CellRef::data() {
    data_t& ret = group->cells[index].data;
    return ret;
}


double CellRef::size(int d) const {
    return gtree.size[d] / (2 << level());
}


double CellRef::origin(int d) const {
    return group->origin[d] + (((index >> d) & 1) * size(d));
}


double CellRef::centre(int d) const {
    return origin(d) + size(d)*0.5;
}

double CellRef::vertex(int v, int d) const {
    return origin(d) + (((v >> d) & 1) * size(d));
}

bool CellRef::containsPoint(double pt[DIM]) const {
    for (int d = 0; d < DIM; d++) {
        if (pt[d] < origin(d) || pt[d] > origin(d)+size(d))
            return false;
    }
    return true;
}


bool CellRef::operator==(const CellRef& rhs) const {
    return group == rhs.group && index == rhs.index;
}

bool CellRef::operator<(const CellRef& rhs) const {
    return group < rhs.group || index < rhs.index;
}

} // namespace tftt


std::ostream& operator<<(std::ostream& os, const tftt::CellRef& cr) {
    if (!cr.isValid())
        os << "{null}";
    else if (cr.isCopyBoundary())
        os << "{copy}";
    else if (cr.isReflectBoundary())
        os << "{reflect}";
    else
        os << "{" << cr.id() << "}";
}

