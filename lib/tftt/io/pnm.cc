
#include <string>
#include <fstream>
#include <ostream>

#include "../config.h"
#include "../cellref.h"
#include "../iter/leaves.h"
#include "../structure/tree.h"

#include "pnm.h"


namespace tftt {
namespace pnm {


void pgm(std::string fname, int imgW, int imgH, fnCell data)
{
    std::ofstream ofs(fname, std::ios::binary);
    pgm(ofs, imgW, imgH, data);
}


void pgm(std::ostream& os, int imgW, int imgH, fnCell data)
{
    unsigned char* bmp = new unsigned char[imgW*imgH];

    int x1,y1,w,h,x,y;
    for (auto& cell : leaves) {
        x1 = int((cell.origin(0)) * (imgW-1) / gtree.size[0]);
        y1 = int((cell.origin(1)) * (imgH-1) / gtree.size[1]);
        w = int(cell.size(0) * (imgW-1) / gtree.size[0])+1;
        h = int(cell.size(1) * (imgH-1) / gtree.size[1])+1;

        for (y = y1; y < y1+h; y++) {
            for (x = x1; x < x1 + w; x++) {
                bmp[y*imgW+x] = (unsigned char)(data(cell)*255.0);
            }
        }
    }

    os << "P5 " << imgW << " " << imgH << " " << 255 << "\n";

    for (y = 0; y < imgH; y++) {
        for (x = 0; x < imgW; x++) {
            os << (unsigned char)(255-bmp[y*imgW+x]);
        }
    }

    delete[] bmp;
}


} // namespace pnm
} // namespace tftt
