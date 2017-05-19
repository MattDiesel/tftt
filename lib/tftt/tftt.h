
#ifndef TFTT_PUBLIC_H
#define TFTT_PUBLIC_H

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

#include "config.h"
#include "treeid.h"

#include "cellref.h"
#include "faceref.h"
#include "vertexref.h"


// bring utilities into tftt namespace
#include "leaves.h"
#include "tfttio.h"
#include "tfttops.h"
#include "adapt.h"
#include "fttcore.h"

#ifndef TFTT_NOMPI
    #include "parallel.h"
#endif


namespace tftt {

// Miscellaneous functions that are mostly deprecated

bool checkAround(cell_t cl, int dist, fnCheckCell check);
bool findAround(cell_t cl, int dist, fnCheckCell check);

double interpFace(cell_t cl, int fc, fnData dt);
double interpChild(cell_t cl, int ch, int forNb, fnData dt);
double interpALEVertex(cell_t cl, int v, fnData dt);


} // namespace tftt


#endif
