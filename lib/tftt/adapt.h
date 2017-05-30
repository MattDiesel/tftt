
#ifndef TFTT_ADAPT_H
#define TFTT_ADAPT_H


#include <set>

#include "config.h"
#include "structure/treecell.h"


namespace tftt {


void twoToOne(cell_t cl);

extern std::set<cell_t, cell_t::less> adaptList;
void adaptBegin();

void adaptAdd(cell_t cr);
bool adaptCommit();

void adaptAddCoarsen(cell_t cr);
bool adaptCommitCoarsen();


void adaptSwBegin();
void adaptSwCommit();
void adaptSwSetCoarsen(cell_t cl);
void adaptSwSetRefine(cell_t cl);
void adaptSwSetHoldRefined(cell_t cl);

void adaptSwPropogateLevel(cell_t cl, int dir, int lvl);
void adaptSwSetFlags(cell_t cl, ADAPTFLAGS af);
ADAPTFLAGS adaptSwGetFlags(cell_t cl);


} // namespace tftt


#endif
