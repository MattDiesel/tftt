
#include <iostream>

#include "tftt/tftt.h"
#include "util/pars.h"

#include "readData.h"
#include "phys.h"
#include "init.h"
#include "compute.h"
#include "solve.h"


using namespace util; // getpars


int t;

void timestep();


int main(int argc, char const *argv[])
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <PAR FILE>" << std::endl;
        return 1;
    }

    try {
        std::cout << "Reading parameters from " << argv[1] << std::endl;
        getpars(argv[1]);
        readData();
    }
    catch (std::exception& e) {
        std::cout << "Error reading parameter file: " << e.what() << std::endl;
    }


    // Init
    physInit();

    initCondition();


    for(t=0; t<=tmax; t++) {
        if(t % tprint == 0) {
            std::cout << "Iteration " << t << ", Time = " << tau*t << '\n';
        }
        timestep();
    }
}


void timestep() {
    double eps;
    tftt::cell_t maxDive;
    // TODO: PLIC for interface tracking

    computeRho();
    computeAdv();
    computeVisc();
    computeIntermVelo();
    computeFlux();
    computeDive();

    if (t % tprint == 0) {
        maxDive = tftt::findmax([](tftt::data_t& dt) { return dt.dive; }, &eps);

        std::cout << "\tBefore: Max div = " << maxDive->dive << " @ " << maxDive << '\n';
        eps = computeCompatability();
        std::cout << "\tCompatibility: " << eps << std::endl;
    }

    computePoisCoef();



    correctVelo();
    correctFlux();
    computeDive();
    if(t % tprint == 0)
    {
        maxDive = tftt::findmax([](tftt::data_t& dt) { return dt.dive; }, &eps);
        std::cout << "\tAfter: Max div = " << maxDive->dive << " @ " << maxDive << '\n';
    }
}
