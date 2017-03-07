
#ifdef TFTT_DEBUG
    #include <stdexcept>
#endif

#include <sstream>
#include <ostream>

#include "tree.h"
#include "treecell.h"
#include "hilbert.h"

#include "cellref.h"


namespace tftt {


CellRef::CellRef()
        : group(nullptr), index(0) {
}

CellRef::CellRef(int flag)
        : group(nullptr), index(flag) {
}


CellRef::CellRef(TreeGroup* gr, int ind)
        : group(gr), index(ind) {
    #ifdef TFTT_DEBUG
        if (!gr) throw std::invalid_argument("Null TreeGroup given.");
        if (ind < 0 || ind > (1 << DIM)) throw std::invalid_argument("Child index is invalid for dimension of problem.");
    #endif
}


bool CellRef::isValid() const {
    return group || index;
}

bool CellRef::isRoot() const {
    return index == -1;
}


bool CellRef::isBoundary() const {
    return (group == gtree.boundGroups || group->isBoundary());
}

int CellRef::boundary() const {
    if (group == gtree.boundGroups) {
        return children()->boundary;
    }
    return group->boundary;
}


CellRef CellRef::parent() const {
    return group->parent;
}

CellRef CellRef::child(int n) const {
    return CellRef(children(), n);
}

CellRef CellRef::childOnFace(int fc, int n) const {
    int start = (fc & 1) << (fc >> 1);
    int step = 1 << (((fc >> 1) + 1) % DIM);

    return child((start + n*step) % (1 << DIM));
}

bool CellRef::hasChildren() const {
    return isRoot() || (group && children());
}

bool CellRef::hasGrandChildren() const {
    if (!hasChildren()) return false;

    for (auto& ch : children()->cells) {
        if (ch.children) return true;
    }
    return false;
}

TreeGroup* CellRef::children() const {
    if (isRoot())
        return gtree.root;
    else
        return group->cells[index].children;
}


// TreeGroup* CellRef::group() const {
//     #ifdef TFTT_DEBUG
//         if (isBoundary()) {
//             throw std::runtime_error("Using boundary group as standard group");
//         }
//     #endif

//     return reinterpret_cast<TreeGroup*>(group);
// }

// TreeGroup* CellRef::bgroup() const {
//     #ifdef TFTT_DEBUG
//         if (!isBoundary()) {
//             throw std::runtime_error("Using standard group as boundary group");
//         }
//     #endif

//     return reinterpret_cast<TreeBoundaryGroup*>(group);
// }

bool nbInParent(int ch, int nb) {
    // Explanation:
    // Reduce to 1d problem, then check bounds.
    
    return ((ch >> (nb >> 1)) ^ nb) & 1;
}


CellRef CellRef::neighbour(int n) const {
    if (nbInParent(index, n))
        return CellRef(group, index ^ (1 << (n >> 1)));

    CellRef cr = group->neighbours[n];
    if (!cr.hasChildren())
        return cr;

    return cr.child(index ^ (1 << (n >> 1)));
}

int CellRef::orientation() const {
    // = f(index, group->orientation)

    return hilbOrient(group->orientation, index);
}


CellRef CellRef::next() const {
    int hch = hilbInvChild(group->orientation, index);
    CellRef ret;

    if (hch == (1<<DIM)-1) { // Last in group
        ret = group->next;

#ifdef TFTT_DEBUG
        if (ret.hasChildren()) {
            throw std::runtime_error("Curve broken, next element has children.");
        }
#endif

        return ret;
    }

    ret = CellRef(group, hilbChild(group->orientation, hch+1));

    while (ret.hasChildren()) {
        ret = ret.child(hilbChild(ret.orientation(), 0));
    }

    return ret;
}


CellRef CellRef::prev() const {
    int hch = hilbInvChild(group->orientation, index);

    if (hch == 0) { // First in group
        return group->prev;
    }

    CellRef ret = CellRef(group, hilbChild(group->orientation, hch-1));

    while (ret.hasChildren()) {
        ret = ret.child(hilbChild(ret.orientation(), (1<<DIM)-1));
    }

    return ret;
}


bool CellRef::isLastInGroup() const {
    return hilbIsLast(group->orientation, index);
}


bool CellRef::isFirstInGroup() const {
    return hilbIsFirst(group->orientation, index);
}


node_t& CellRef::rank() {
    return group->cells[index].rank;
}

node_t const& CellRef::rank() const {
    return group->cells[index].rank;
}


ident_t CellRef::id() const {
    if (isRoot()) return 0;
    return group->id.id | index;
}

int CellRef::level() const {
    if (isRoot()) return -1; // Todo: Needed? 
    return group->id.level();
}


data_t& CellRef::data() {
#ifdef TFTT_DEBUG
    if (hasChildren()) {
        throw std::runtime_error("Attempt to access data from non-leaf");
    }
#endif

    data_t& ret = group->cells[index].data;
    return ret;
}

data_t const& CellRef::data() const {
    if (hasChildren()) {
        throw std::runtime_error("Attempt to access data from non-leaf");
    }

    data_t& ret = group->cells[index].data;
    return ret;
}


double CellRef::avrChildren(fnData dt) const {
    double ret = 0.0;
    for (int ch = 0; ch < 1<<DIM; ch++) {
        if (child(ch).hasChildren())
            ret += child(ch).avrChildren(dt);
        else
            ret += dt(child(ch).data());
    }
    return ret / (1<<DIM);
}


double CellRef::ngbVal(int nb, fnData dt, double* ifBoundary) const {
    cell_t ngb = neighbour(nb);

    // if (ngb.isBoundary()) {
    //     if (ifBoundary)
    //         return *ifBoundary;
    //     else
    //         throw std::runtime_error("Not Implemented: Automatic boundary condition");
    // }
    // else
    if (ngb.hasChildren()) {
        // Neighbour more refined, average children
        return ngb.avrChildren(dt);
    }
    else if (ngb.level() < level()) {
        // Neighbour is less refined
        return interpChild(ngb, index ^ (1 << (nb >> 1)), nb ^ 1, dt);
    }
    else {
        // Neighbour at same level
        return dt(ngb.data());
    }
}


// facedata_t& CellRef::facedata(int dir) {
//     facedata_t& ret = group->cells[index].facedata[dir];
//     return ret;
// }

// facedata_t const& CellRef::facedata(int dir) const {
//     facedata_t& ret = group->cells[index].facedata[dir];
//     return ret;
// }


face_t CellRef::face(int dir) const {
    return group->cells[index].faces[dir];
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


void CellRef::sizes(double ret[DIM]) const {
    int lvl = level();

    for (int d = 0; d < DIM; d++) {
        ret[d] = gtree.size[d] / (2 << lvl);
    }
}

void CellRef::origins(double ret[DIM]) const {
    int lvl = level();

    for (int d = 0; d < DIM; d++) {
        ret[d] = group->origin[d]
            + (((index >> d) & 1) * (gtree.size[d] / (2 << lvl)));
    }
}

void CellRef::vertices(double ret[1<<DIM][DIM]) const {
    int lvl = level();
    double siz;
    double org;

    for (int d = 0; d < DIM; d++) {
        siz = gtree.size[d] / (2 << lvl);
        org = group->origin[d] + (((index >> d) & 1) * siz);

        for (int v = 0; v < 1 << DIM; v++) {
            ret[v][d] = org + (((v >> d) & 1) * siz);
        }
    }
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


bool CellRef::operator!=(const CellRef& rhs) const {
    return group != rhs.group || index != rhs.index;
}


bool CellRef::operator<(const CellRef& rhs) const {
    return group < rhs.group || index < rhs.index;
}


data_t* CellRef::operator->() {
    if (hasChildren()) {
        throw std::runtime_error("Attempt to access data from non-leaf");
    }
    
    return &group->cells[index].data;
}

data_t const* CellRef::operator->() const {
    if (hasChildren()) {
        throw std::runtime_error("Attempt to access data from non-leaf");
    }

    return &group->cells[index].data;
}


std::string CellRef::path() const {
    std::ostringstream oss;
    oss << *this;
    return oss.str();
}


} // namespace tftt


std::ostream& operator<<(std::ostream& os, const tftt::CellRef& cr) {
    if (!cr.isValid())
        os << "{null}";
    else if (cr.isBoundary())
        os << "{boundary(" << cr.boundary() << ") " << cr.id() << "}";
    else
        os << "{" << cr.id() << "}";
}

