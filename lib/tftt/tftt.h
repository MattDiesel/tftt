

#ifndef TFTT_PUBLIC_H
#define TFTT_PUBLIC_H


#include <cstdint>
#include <string>
#include <ostream>
#include <set>


#ifndef DIM
#define DIM 2
#endif


// Constants for children, vertices, neighbours

// Neighbours. P = positive dir, N = negative
enum {
    NB_XN = 0,
    NB_XP,
    NB_YN,
    NB_YP,
    NB_ZN,
    NB_ZP
};

// Vertices (2d)
enum {
    VT_BL = 0,
    VT_BR = 1,
    VT_TL = 2,
    VT_TR = 3
}

// Children (2d)
enum {
    CH_BL = 0,
    CH_BR = 1,
    CH_TL = 2,
    CH_TR = 3
}



struct rt_data {
    double cc;
    // double vof1;
    // double vof2;
    // double vof3;
    double P;
    double Pstar;
    double rhs;
    double cenCoef;
    double res;

    // double vof;
    double rho;
    double V[DIM];
    double v[DIM];
    double F1[DIM];
    double F2[DIM];
    double D1[DIM];
    double D2[DIM];

    double dive;
};

struct rt_facedata {
    double poisCoef;
};

namespace tftt {
	template<typename T> struct TreeId;

    typedef rt_data data_t;
    typedef rt_facedata facedata_t;
	typedef TreeId<uint64_t> ident_t;


    // Functions to be supplied by user
    typedef double (*fnDataNorm)(data_t& d, int max);
    typedef double (*fnData)(data_t& d);

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
cell_t findmax(fnData dt);

double interpChild(cell_t cl, int ch, int forNb, fnData dt);
double interpALEVertex(cell_t cl, int v, fnData dt);

cell_t atPos(double pos[DIM]);
cell_t atVertex(int v);

void refine(CellRef cl);
void coarsen(CellRef cl);
void twoToOne(CellRef cl);


extern std::set<cell_t> adaptList;
void adaptBegin();
void adaptAdd(CellRef cr);
bool adaptCommit();

void adaptAddCoarsen(CellRef cr);
bool adaptCommitCoarsen();


} // namespace tftt


// bring utilities into tftt namespace
#include "leaves.h"
#include "tfttio.h"


#endif
