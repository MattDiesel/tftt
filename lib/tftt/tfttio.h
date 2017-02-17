
#ifndef TFTT_IO_H
#define TFTT_IO_H


#include <string>
#include <ostream>

#include "tftt.h"


namespace tftt {

void drawCell(std::ostream& os, cell_t const& c);

void drawMesh(std::string fname);
void drawMesh(std::ostream& os);


void drawCurve(std::string fname);
void drawCurve(std::ostream& os);

void drawBoundaries(std::string fname);
void drawBoundary(std::ostream& os, int b);

void drawMatrix(std::string fname, int imgW, int imgH, fnDataNorm dataNorm);
void drawMatrix(std::ostream& os, int imgW, int imgH, fnDataNorm dataNorm);


void saveTree(std::string fname);
void saveTree(std::ostream& os);

void loadTree(std::string fname);


void splitToDisk(std::string fnameFmt);


} // namespace tftt

#endif
