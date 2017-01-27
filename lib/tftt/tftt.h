

#ifndef TFTT_PUBLIC_H
#define TFTT_PUBLIC_H


#include <cstdint>
#include <string>
#include <ostream>


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

cell_t find(ident_t id);
void refine(CellRef cl);


void drawMesh(std::string fname);
void drawMesh(std::ostream& os);



} // namespace tftt


#endif
