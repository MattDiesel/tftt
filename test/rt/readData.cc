
#include <stdexcept>

#include "util/pars.h"

#include "phys.h"

#include "readData.h"


using namespace util; // tfetch


// #define Q(x) #x
// #define QUOTE(x) Q(x)
// #define checkAndSet(x) if(!tfetch(QUOTE(x), x)) { \
//         throw std::out_of_range("missing " QUOTE(x)); }


double tau, L[2], delta[2];
int N[2], tmax, tprint, tvelo;
int wingSize[2];
double headPos[2];


void readData()
{
    if(!tfetch("tau", tau)) 
        throw std::out_of_range("missing tau");

    if(!tfetch("Lx", L[0]))
        throw std::out_of_range("missing Lx");
    if(!tfetch("Ly", L[1]))
        throw std::out_of_range("missing Ly");
    if(!tfetch("nx", N[0]))
        throw std::out_of_range("missing Nx");
    if(!tfetch("ny", N[1]))
        throw std::out_of_range("missing Ny");

    delta[0] = L[0]/N[0];
    delta[1] = L[1]/N[1];

    // Boundary
    wingSize[0] = N[0]/10;
    wingSize[1] = N[1]/5;
    tfetch("wingHeight", wingSize[0]);
    tfetch("wingWidth", wingSize[1]);
    headPos[0] = (L[0]/4);
    tfetch("xHead", headPos[0]);
    headPos[1] = L[1]/2;
    tfetch("yHead", headPos[1]);

    if(!tfetch("tmax", tmax)) 
        throw std::out_of_range("missing tmax");
    if(!tfetch("tprint", tprint)) 
        throw std::out_of_range("missing tprint");

    tvelo = tprint;
    tfetch("tvelo", tvelo);

    physReadData();
}
