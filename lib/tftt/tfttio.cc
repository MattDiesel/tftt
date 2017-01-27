
#include <string>
#include <fstream>
#include <ostream>

#include "tftt.h"
#include "tree.h"
#include "gray.h"

#include "tfttio.h"


namespace tftt {


void drawMesh(std::string fname) {
    std::ofstream ofs(fname);
    drawMesh(ofs);
}

void drawMesh(std::ostream& os) {
    for (auto& c : leaves) {
        for (int v = 0; v < 1<<DIM; v++) {
            os << c.vertex(utils::toGray(v), 0);
            for (int d = 1; d < DIM; d++) {
                os << " " << c.vertex(utils::toGray(v), d);
            }
            os << "\n";
        }
        os << "\n";
    }
}


void drawMatrix(std::string fname, int imgW, int imgH) {
    std::ofstream ofs(fname, std::ios::binary);
    drawMatrix(ofs, imgW, imgH);
}

void drawMatrix(std::ostream& os, int imgW, int imgH) {
    unsigned char* bmp = new unsigned char[imgW*imgH];

    int x1,y1,w,h,x,y;
    for (auto& cell : leaves) {
        x1 = int((cell.origin(0)) * (imgW-1) / gtree.size[0]);
        y1 = int((cell.origin(1)) * (imgH-1) / gtree.size[1]);
        w = int(cell.size(0) * (imgW-1) / gtree.size[0])+1;
        h = int(cell.size(1) * (imgH-1) / gtree.size[1])+1;

        for (y = y1; y < y1+h; y++) {
            for (x = x1; x < x1 + w; x++) {
                bmp[y*imgW+x] = (unsigned char)(cell.data() * 255.0);
            }
        }
    }

    os << "P5 " << imgW << " " << imgH << " " << 255 << "\n";

    for (y = 0; y < imgH; y++) {
        for (x = 0; x < imgW; x++) {
            os << (unsigned char)(255-bmp[y*imgW+x]);
        }
    }
}


} // namespace tftt
