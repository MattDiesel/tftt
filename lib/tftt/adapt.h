
#ifndef TFTT_ADAPT_H
#define TFTT_ADAPT_H


#include <set>

#include "config.h"
#include "structure/treecell.h"
#include "cellref.h"


namespace tftt {


void twoToOne(CellRef cl);

extern std::set<CellRef, CellRef::less> adaptList;
void adaptBegin();

void adaptAdd(CellRef cr);
bool adaptCommit();

void adaptAddCoarsen(CellRef cr);
bool adaptCommitCoarsen();


void adaptSwBegin();
void adaptSwCommit();
void adaptSwSetCoarsen(CellRef cl);
void adaptSwSetRefine(CellRef cl);
void adaptSwSetHoldRefined(CellRef cl);

void adaptSwPropogateLevel(CellRef cl, int dir, int lvl);
void adaptSwSetFlags(CellRef cl, ADAPTFLAGS af);
ADAPTFLAGS adaptSwGetFlags(CellRef cl);


} // namespace tftt


#endif
