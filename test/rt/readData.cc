
#include <stdexcept>

#include "pars.h"
#include "phys.h"

#include "readData.h"

// #define Q(x) #x
// #define QUOTE(x) Q(x)
// #define dcheckAndSet(x) if(!dfetch(QUOTE(x), x)) { \
//         throw std::out_of_range("missing " QUOTE(x)); }
// #define icheckAndSet(x) if(!ifetch(QUOTE(x), x)) { \
//         throw std::out_of_range("missing " QUOTE(x)); }


double tau, L[2], delta[2];
int N[2], tmax, tprint, tvelo;
int wingSize[2];
double headPos[2];


void readData()
{
    if(!dfetch("tau", tau)) 
        throw std::out_of_range("missing tau");

    if(!dfetch("Lx", L[0]))
        throw std::out_of_range("missing Lx");
    if(!dfetch("Ly", L[1]))
        throw std::out_of_range("missing Ly");
    if(!ifetch("nx", N[0]))
        throw std::out_of_range("missing Nx");
    if(!ifetch("ny", N[1]))
        throw std::out_of_range("missing Ny");

    delta[0] = L[0]/N[0];
    delta[1] = L[1]/N[1];

    // Boundary
    wingSize[0] = N[0]/10;
    wingSize[1] = N[1]/5;
    ifetch("wingHeight", wingSize[0]);
    ifetch("wingWidth", wingSize[1]);
    headPos[0] = (L[0]/4);
    dfetch("xHead", headPos[0]);
    headPos[1] = L[1]/2;
    dfetch("yHead", headPos[1]);

    if(!ifetch("tmax", tmax)) 
        throw std::out_of_range("missing tmax");
    if(!ifetch("tprint", tprint)) 
        throw std::out_of_range("missing tprint");

    tvelo = tprint;
    ifetch("tvelo", tvelo);

    physReadData();
}
