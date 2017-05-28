
#ifndef TFTT_STRUCTURE_GROUPMEMORY_H
#define TFTT_STRUCTURE_GROUPMEMORY_H


namespace tftt {

struct TreeGroup;
struct TreeCell;

namespace group {

template<class ...Us>
TreeGroup* create(Us... args)
{
    return new TreeGroup(args...);
}

void free(TreeGroup* grp)
{
    delete grp;
}

int getCellIndex(TreeCell* tc)
{
    return tc->index;
}

TreeGroup* getCellGroup(TreeCell* tc)
{
    return tc->group;
}


} // namespace group
} // namespace tftt


#endif
