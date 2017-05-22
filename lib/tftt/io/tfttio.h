
#ifndef TFTT_IO_H
#define TFTT_IO_H


#include <string>
#include <ostream>

#include "../config.h"
#include "../cellref.h"


namespace tftt {


void drawCell(std::ostream& os, cell_t const& c);

void drawMesh(std::string fname);
void drawMesh(std::ostream& os);


void drawPartialMesh(std::string fname);
void drawPartialMesh(std::string fname, cell_t from, cell_t to);
void drawPartialMesh(std::ostream& os, cell_t from, cell_t to);


void drawGhosts(std::string fname, int n);
void drawGhosts(std::ostream& os, int n);

void drawBorder(std::string fname, int n);
void drawBorder(std::ostream& os, int n);


void drawCurve(std::string fname);
void drawCurve(std::ostream& os);

void drawPartialCurve(std::string fname);
void drawPartialCurve(std::ostream& os);

void drawBoundaries(std::string fname);
void drawBoundary(std::ostream& os, int b);

void drawMatrix(std::string fname, int imgW, int imgH, fnDataNorm dataNorm);
void drawMatrix(std::ostream& os, int imgW, int imgH, fnDataNorm dataNorm);

void plotMatrix(std::string fname, fnData dt);
void plotMatrix(std::ostream& os, fnData dt);


void drawPoissonNeighbourhood(std::string fname, cell_t cl);
void drawPoissonNeighbourhood(std::ostream& os, cell_t cl);


} // namespace tftt

#endif
