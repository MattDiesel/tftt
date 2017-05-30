
#include <string>
#include <fstream>
#include <ostream>

#include "../config.h"
#include "../cellref.h"
#include "../fttcore.h"
#include "../iter/activecurve.h"
#include "../structure/tree.h"

#include "plot3d.h"


namespace tftt {
namespace plot3d {


void scatter(std::string fname, fnCell dt) {
    std::ofstream ofs(fname);
    scatter(ofs, dt);
}


void scatter(std::ostream& os, fnCell dt) {
    for (auto& cl : activecurve) {
        os << cl.centre(0) << " " << cl.centre(1) << " " << dt(cl) << "\n";
    }
}


void sample(std::string fname, int resX, int resY, fnCell dt)
{
    std::ofstream ofs(fname);
    sample(ofs, resX, resY, dt);
}


void sample(std::ostream& os, int resX, int resY, fnCell dt)
{
    double dX = gtree.size[0] / resX;
    double dY = gtree.size[1] / resY;

    cell_t cl;
    double p[2];
    for (int y = 0; y < resY; y++) {
        for (int x = 0; x < resX; x++) {
            p[0] = x*dX+dX*0.5;
            p[1] = y*dY+dY*0.5;
            cl = tftt::atPos(p);

            os << p[0] << " " << p[1] << " " << dt(cl) << "\n";
        }
        os << "\n";
    }
}


} // namespace plot
} // namespace tftt
