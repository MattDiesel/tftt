
#include <stdexcept>
#include <string>
#include <sstream>
#include <ostream>

#include "structure/tree.h"
#include "structure/treecell.h"
#include "structure/treegroup.h"
#include "structure/groupmemory.h"
#include "hilbert.h"
#include "iter/all.h"

#include "cellref.h"


namespace tftt {


inline TreeGroup* cell_t::group() const
{
    return group::getCellGroup(_cell);
}

inline int cell_t::index() const
{
    return group::getCellIndex(_cell);
}

inline TreeCell* cell_t::treecell() const
{
    return _cell;
}

inline void cell_t::stepChild()
{
    _cell++;
}


inline CellRef<TreeCell>::CellRef()
    : _cell(nullptr)
{
}

inline CellRef<TreeCell>::CellRef(int flag)
    : _cell(nullptr)
{
    if (flag == -1)
        _cell = &rootCell;
}


inline CellRef<TreeCell>::CellRef(TreeGroup* gr, int ind)
{
    _cell = &gr->cells[ind];
    // #ifdef TFTT_DEBUG
    // if (!gr) throw std::invalid_argument("Null TreeGroup given.");
    // if ((gr != gtree.root && ind < 0)
    //         || ind > (1 << DIM)) throw std::invalid_argument("Child index is invalid for dimension of problem.");
    // #endif
}


inline bool cell_t::isValid() const
{
    return treecell();
}

inline bool cell_t::isRoot() const
{
    return treecell() == &rootCell;
}


inline bool cell_t::isBoundary() const
{
    return (group() == gtree.boundGroups || group()->isBoundary());
}

inline int cell_t::boundary() const
{
    if (group() == gtree.boundGroups) {
        return index();
    }
    return group()->boundary();
}


inline cell_t cell_t::parent() const
{
    return group()->parent;
}

inline cell_t cell_t::child(int n) const
{
    return cell_t(children(), n);
}

inline cell_t cell_t::childOnFace(int fc, int n) const
{
    return child(childIndexOnFace(fc, n));
}

inline int cell_t::childIndexOnFace(int fc, int n)
{
    int start = (fc & 1) << (fc >> 1);
    int step = 1 << (((fc >> 1) + 1) % DIM);

    return ((start + n*step) % (1 << DIM));
}

inline bool cell_t::hasChildren() const
{
    return isRoot() || (group() && children());
}

inline bool cell_t::hasGrandChildren() const
{
    if (!hasChildren()) return false;

    for (auto& ch : children()->cells) {
        if (ch.children) return true;
    }
    return false;
}

inline TreeGroup* cell_t::children() const
{
    if (isRoot())
        return gtree.root;
    else
        return treecell()->children;
}


// TreeGroup* cell_t::group()() const {
//     #ifdef TFTT_DEBUG
//         if (isBoundary()) {
//             throw std::runtime_error("Using boundary group() as standard group()");
//         }
//     #endif

//     return reinterpret_cast<TreeGroup*>(group());
// }

// TreeGroup* cell_t::bgroup()() const {
//     #ifdef TFTT_DEBUG
//         if (!isBoundary()) {
//             throw std::runtime_error("Using standard group() as boundary group()");
//         }
//     #endif

//     return reinterpret_cast<TreeBoundaryGroup*>(group());
// }

inline bool nbInParent(int ch, int nb)
{
    // Explanation:
    // Reduce to 1d problem, then check bounds.

    return ((ch >> (nb >> 1)) ^ nb) & 1;
}


inline cell_t cell_t::neighbour(int n) const
{
    int ind = index();

    if (nbInParent(ind, n)) {
        // return cell_t(group(), index() ^ (1 << (n >> 1)));

        int diff = (ind ^ (1 << (n >> 1))) - ind;
        cell_t ret(*this);

        ret._cell += diff;
        return ret;
    }

    cell_t cr = group()->neighbours[n];
    if (!cr.hasChildren())
        return cr;

    return cr.child(index() ^ (1 << (n >> 1)));
}


inline cell_t cell_t::diagonal(int n) const
{
    int ch = ~index() & (2*DIM-1);
    cell_t cr = *this;

    if (ch == n) {
        // In parent
        // cr = cell_t(group(), ch);
        cr = cell_t(*this);
        cr._cell += (n - index());

        if (cr.hasChildren())
            cr = cr.child(~n & (2*DIM-1));

        return cr;
    }

    int x;
    for (int i = 0; i < DIM; i++) {
        x = (i << 1) | ((n >> i) & 1);
        cr = cr.neighbour(x);

        if (cr.hasChildren())
            cr = cr.child(~n & (2*DIM-1));
    }

    return cr;
}


inline tagNeighbours cell_t::neighbours() const
{
    return tagNeighbours(*this);
}


inline tagPoissonNeighbours cell_t::poissonNeighbourhood() const
{
    return tagPoissonNeighbours(*this);
}


inline int cell_t::orientation() const
{
    // = f(index(), group()->orientation)

    return hilbOrient(group()->orientation, index());
}


inline cell_t cell_t::next() const
{
    int hch = hilbInvChild(group()->orientation, index());
    cell_t ret;

    if (hch == (1<<DIM)-1) { // Last in group()
        ret = group()->next;

        #ifdef TFTT_DEBUG
        if (ret.isValid() && ret.hasChildren()) {
            throw std::runtime_error("Curve broken, next element has children.");
        }
        #endif

        return ret;
    }

    ret = cell_t(group(), hilbChild(group()->orientation, hch+1));

    while (ret.hasChildren()) {
        ret = ret.child(hilbChild(ret.orientation(), 0));
    }

    return ret;
}


inline cell_t cell_t::prev() const
{
    int hch = hilbInvChild(group()->orientation, index());

    if (hch == 0) { // First in group()
        return group()->prev;
    }

    cell_t ret = cell_t(group(), hilbChild(group()->orientation, hch-1));

    while (ret.hasChildren()) {
        ret = ret.child(hilbChild(ret.orientation(), (1<<DIM)-1));
    }

    return ret;
}


inline bool cell_t::isLastInGroup() const
{
    return hilbIsLast(group()->orientation, index());
}


inline bool cell_t::isFirstInGroup() const
{
    return hilbIsFirst(group()->orientation, index());
}


inline node_t& cell_t::rank()
{
    return treecell()->rank;
}

inline node_t const& cell_t::rank() const
{
    return treecell()->rank;
}


inline ident_t cell_t::id() const
{
    if (isRoot()) return 0;
    return group()->id.id | index();
}

inline int cell_t::level() const
{
    if (isRoot()) return -1; // Todo: Needed?
    return group()->id.level();
}


inline data_t& cell_t::data()
{
    #ifdef TFTT_DEBUG
    if (hasChildren()) {
        throw std::runtime_error("Attempt to access data from non-leaf");
    }
    #endif

    data_t& ret = treecell()->data;
    return ret;
}

inline data_t const& cell_t::data() const
{
    if (hasChildren()) {
        throw std::runtime_error("Attempt to access data from non-leaf");
    }

    data_t& ret = treecell()->data;
    return ret;
}


inline double cell_t::avrChildren(fnConstCell dt) const
{
    double ret = 0.0;
    if (!hasChildren()) return dt(*this);

    for (int ch = 0; ch < 1<<DIM; ch++) {
        if (child(ch).hasChildren())
            ret += child(ch).avrChildren(dt);
        else
            ret += dt(child(ch));
    }
    return ret / (1<<DIM);
}


inline double cell_t::ngbVal(int nb, fnData dt, double* ifBoundary) const
{
    throw std::runtime_error("Deprecated.");
    return 0.0;

    // cell_t ngb = neighbour(nb);

    // // if (ngb.isBoundary()) {
    // //     if (ifBoundary)
    // //         return *ifBoundary;
    // //     else
    // //         throw std::runtime_error("Not Implemented: Automatic boundary condition");
    // // }
    // // else
    // if (ngb.hasChildren()) {
    //     // Neighbour more refined, average children
    //     return ngb.avrChildren(dt);
    // }
    // else if (ngb.level() < level()) {
    //     // Neighbour is less refined
    //     return interpChild(ngb, index() ^ (1 << (nb >> 1)), nb ^ 1, dt);
    // }
    // else {
    //     // Neighbour at same level
    //     return dt(ngb.data());
    // }
}

#ifdef TFTT_FACES

inline // facedata_t& cell_t::facedata(int dir) {
//     facedata_t& ret = treecell()facedata[dir];
//     return ret;
// }

inline // facedata_t const& cell_t::facedata(int dir) const {
//     facedata_t& ret = treecell()facedata[dir];
//     return ret;
// }


inline face_t cell_t::face(int dir) const
{
    return treecell()->faces[dir];
}

#endif


#ifdef TFTT_VERTICES

inline vertex_t cell_t::vertex(int v) const
{
    return treecell()->vertices[v];
}

#endif


inline double cell_t::size(int d) const
{
    if (isRoot()) return gtree.size[d];
    return gtree.size[d] / (2 << level());
}


inline double cell_t::origin(int d) const
{
    if (isRoot()) return 0.0;
    return group()->origin[d] + (((index() >> d) & 1) * size(d));
}


inline double cell_t::centre(int d) const
{
    return origin(d) + size(d)*0.5;
}

inline double cell_t::vertexPoint(int v, int d) const
{
    return origin(d) + (((v >> d) & 1) * size(d));
}


inline void cell_t::sizes(double ret[DIM]) const
{
    int lvl = level();

    for (int d = 0; d < DIM; d++) {
        ret[d] = gtree.size[d] / (2 << lvl);
    }
}

inline void cell_t::origins(double ret[DIM]) const
{
    int lvl = level();

    for (int d = 0; d < DIM; d++) {
        ret[d] = group()->origin[d]
                 + (((index() >> d) & 1) * (gtree.size[d] / (2 << lvl)));
    }
}

inline void cell_t::vertices(double ret[1<<DIM][DIM]) const
{
    int lvl = level();
    double siz;
    double org;

    for (int d = 0; d < DIM; d++) {
        siz = gtree.size[d] / (2 << lvl);
        org = group()->origin[d] + (((index() >> d) & 1) * siz);

        for (int v = 0; v < 1 << DIM; v++) {
            ret[v][d] = org + (((v >> d) & 1) * siz);
        }
    }
}



inline bool cell_t::containsPoint(double pt[DIM]) const
{
    for (int d = 0; d < DIM; d++) {
        if (pt[d] < origin(d) || pt[d] > origin(d)+size(d))
            return false;
    }
    return true;
}


inline bool cell_t::operator==(const cell_t& rhs) const
{
    return treecell() == rhs.treecell();
}


inline bool cell_t::operator!=(const cell_t& rhs) const
{
    return treecell() != rhs.treecell();
}


inline bool cell_t::operator<(const cell_t& rhs) const
{
    return treecell() < rhs.treecell();
}


inline data_t* cell_t::operator->()
{
    if (hasChildren()) {
        throw std::runtime_error("Attempt to access data from non-leaf");
    }

    return &treecell()->data;
}

inline data_t const* cell_t::operator->() const
{
    if (hasChildren()) {
        throw std::runtime_error("Attempt to access data from non-leaf");
    }

    return &treecell()->data;
}


inline std::string cell_t::path() const
{
    std::ostringstream oss;

    if (!isValid())
        oss << "{null}";
    else if (isBoundary())
        oss << "{boundary(" << boundary() << ") " << id() << "}";
    else {
        // oss << "{" << cr.id() << "}";
        oss << "{" << centre(0) << "," << centre(1) << "}";
    }

    return oss.str();
}


} // namespace tftt
