

#ifndef TFTT_PUBLIC_H
#define TFTT_PUBLIC_H


#include <cstdint>
#include <string>
#include <ostream>
#include <set>


#ifndef DIM
#define DIM 2
#endif

struct deltawing_data {
    double cc;
    double vof1;
    double vof2;
    double vof3;
    double rho;
    double u;
    double v;
    double p;
    double U;
    double F11;
    double F21;
    double D11;
    double D21;
    double V;
    double F12;
    double F22;
    double D12;
    double D22;
    double dive;
    double res;
    double cenCoef;
    double xCoef;
    double yCoef;
};

namespace tftt {
	template<typename T> struct TreeId;

	typedef deltawing_data data_t;
	typedef TreeId<uint64_t> ident_t;


    // Functions to be supplied by user
    typedef double (*fnDataNorm)(data_t& d, int max);

    extern fnDataNorm dataNorm;
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

cell_t atPos(double pos[DIM]);
cell_t atVertex(int v);

void refine(CellRef cl);
void coarsen(CellRef cl);
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
