
#include <stdexcept>
#include <iostream>

#include "tftt/tftt.h"
#include "util/pars.h"

#include "readData.h"

#include "phys.h"


using namespace util; // tfetch


// Globals
double rho1, rho2, sigma, mu1, mu2, g[2], omega, MinEps;
int npre, npost;


void physReadData() {
    if(!tfetch("rho1", rho1)) 
        throw std::out_of_range("missing rho1");
    if(!tfetch("rho2", rho2)) 
        throw std::out_of_range("missing rho2");
    if(!tfetch("mu1", mu1)) 
        throw std::out_of_range("missing mu1");
    if(!tfetch("mu2", mu2)) 
        throw std::out_of_range("missing mu2");

    // Gravity defaults to 0.0
    g[0] = 0.0;
    g[1] = 9.81;
    tfetch("gravity", g[1]);

    omega = 1.3;
    tfetch("omega", omega);

    MinEps = 1.0e-6;
    tfetch("MinEps", MinEps);
    
    tfetch("npre", npre);
    tfetch("npost", npost);
}


void physInit(void)
{
    std::cout << "omega = " << omega
            << "   MinEps = " << MinEps
            << "   npre = " << npre
            << "   npost = " << npost << std::endl;

    int ng = 0;


    tftt::init(L[0], L[1]);

    // Approximate static mesh
    // int nnx = N[0];
    // int nny = N[1];
    // while (nnx % 2 == 0 && nny % 2 == 0 & nnx > 1 & nny > 1) {
    //     nnx /= 2;
    //     nny /= 2;
    //     ng++;

    //     for (auto& cl : tftt::leaves) {
    //         tftt::refine(cl);
    //     }
    // }

    std::cout << "ng = " << ng << std::endl;

    tftt::drawMesh("initMesh.dat");
}
