

#ifndef TFTT_PUBLIC_H
#define TFTT_PUBLIC_H


#include <cstdint>
#include <string>
#include <ostream>
#include <set>


#ifndef DIM
#define DIM 2
#endif

namespace tftt {
	template<typename T> struct TreeId;

	typedef double data_t;
	typedef TreeId<uint64_t> ident_t;
}


#include "cellref.h" // Only public interface
#include "treeid.h"


namespace tftt {

struct CellRef;
typedef CellRef cell_t;


//! Constructs the top level tree
void init(double w, double h);
void reset();

cell_t find(ident_t id);

void refine(CellRef cl);
void twoToOne(CellRef cl);


extern std::set<cell_t> adaptList;
void adaptBegin();
void adaptAdd(CellRef cr);
bool adaptCommit();


} // namespace tftt


// bring utilities into tftt namespace
#include "leaves.h"
#include "tfttio.h"


#endif
