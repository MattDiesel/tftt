
#ifndef TFTT_STRUCTURE_GROUPMEMORY_H
#define TFTT_STRUCTURE_GROUPMEMORY_H


#include <cstddef>


namespace tftt {

struct TreeGroup;
struct TreeCell;

namespace group {

template<class ...Us>
inline TreeGroup* create(Us... args)
{
    return new TreeGroup(args...);
}

inline void free(TreeGroup* grp)
{
    delete grp;
}

inline int getCellIndex(TreeCell* tc)
{
    return tc->index;
}

inline TreeGroup* getCellGroup(TreeCell* tc)
{
    return reinterpret_cast<TreeGroup*>(
               reinterpret_cast<char*>(tc) - (tc->index*sizeof(TreeCell)) - offsetof(TreeGroup, cells));
}


} // namespace group
} // namespace tftt


#endif
