
#include <stdexcept>
#include <iostream>

#include "tftt/tftt.h"

#include "pars.h"
#include "readData.h"

#include "phys.h"


// Globals
double rho1, rho2, sigma, mu1, mu2, g[2], omega, MinEps;
int npre, npost;


void physReadData() {
    if(!dfetch("rho1", rho1)) 
        throw std::out_of_range("missing rho1");
    if(!dfetch("rho2", rho2)) 
        throw std::out_of_range("missing rho2");
    if(!dfetch("mu1", mu1)) 
        throw std::out_of_range("missing mu1");
    if(!dfetch("mu2", mu2)) 
        throw std::out_of_range("missing mu2");

    // Gravity defaults to 0.0
    g[0] = 0.0;
    g[1] = 9.81;
    dfetch("gravity", g[1]);

    omega = 1.3;
    dfetch("omega", omega);

    MinEps = 1.0e-6;
    dfetch("MinEps", MinEps);
    
    ifetch("npre", npre);
    ifetch("npost", npost);
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
