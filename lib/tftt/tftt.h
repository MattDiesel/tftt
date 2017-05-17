

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
};

// Children (2d)
enum {
    CH_BL = 0,
    CH_BR = 1,
    CH_TL = 2,
    CH_TR = 3
};



struct rt_data {
    double cc;
    // double vof1;
    // double vof2;
    // double vof3;
    double P;
    // double Pstar;
    double rhs;
    double res;

    // double vof;
    // double rho;

    // double dive;
};

struct rt_facedata {
    double poisCoef;
};

struct rt_vertexdata {
    int placeholder;
    // double V[DIM];
    // double v[DIM];
    // double F1[DIM];
    // double F2[DIM];
    // double D1[DIM];
    // double D2[DIM];
};

namespace tftt {

template<typename T> struct TreeId;

typedef rt_data data_t;
typedef rt_facedata facedata_t;
typedef rt_vertexdata vertexdata_t;
typedef TreeId<uint64_t> ident_t;
typedef int8_t node_t;

// Functions to be supplied by user
typedef double (*fnDataNorm)(data_t& d, int max);
typedef double (*fnData)(data_t& d);
typedef double& (*fnDataRef)(data_t& d);

extern struct TFTTOPTIONS {
    int ghostsFlag; // 0 - Minimum, 1 - Children of neighbour groups
    int two2oneFlag; // 0 - None, 1 - strict, 2 - incl. corners, 3 - 3-2-1
    bool isNeuman;
    int maxDepth;
} options;

}


#include "cellref.h" // Todo:Only public interface
#include "treeid.h"
#include "faceref.h"
#include "vertexref.h"


namespace tftt {

struct CellRef;
typedef CellRef cell_t;
typedef FaceRef face_t;
typedef VertexRef vertex_t;

typedef bool (*fnCheckCell)(cell_t& cl);
typedef double (*fnCell)(cell_t& cl);

//! Constructs the top level tree
void init(double w, double h);
void reset();

bool checkAround(cell_t cl, int dist, fnCheckCell check);
bool findAround(cell_t cl, int dist, fnCheckCell check);

double interpFace(cell_t cl, int fc, fnData dt);
double interpChild(cell_t cl, int ch, int forNb, fnData dt);
double interpALEVertex(cell_t cl, int v, fnData dt);

//! Distribute cells across processors
void distribute(int n);

cell_t atPos(double pos[DIM]);
cell_t atVertex(int v);

void refine(CellRef cl);
void coarsen(CellRef cl);
void twoToOne(CellRef cl);


struct crless {
    bool operator()(const CellRef& a, const CellRef& b) {
        if (a.level() < b.level())
            return true;
        else if (a.level() > b.level())
            return false;
        else if (a.group < b.group)
            return true;
        else if (a.group > b.group)
            return false;
        else if (a.index < b.index)
            return true;
        else
            return false;
    }
};

struct crparless {
    bool operator()(const CellRef& a, const CellRef& b) {
        return a.id().id < b.id().id;
    }
};

extern std::set<CellRef, crless> adaptList;
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


} // namespace tftt


// bring utilities into tftt namespace
#include "leaves.h"
#include "tfttio.h"
#include "tfttops.h"

#ifndef TFTT_NOMPI
    #include "parallel.h"
#endif


#endif
