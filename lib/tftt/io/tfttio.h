
#ifndef TFTT_IO_H
#define TFTT_IO_H


#include <string>
#include <ostream>

#include "../config.h"
#include "../cellref.h"


namespace tftt {


//! Writes the coordinates of each cell vertex to a stream
//! \details When opened in gnuplot, this will outline the cell.
void drawCell(std::ostream& os, cell_t const& c);

//! Draws the leaves in the tree to a file
void drawMesh(std::string fname);

//! \copydoc drawMesh()
void drawMesh(std::ostream& os);


#if DIM == 2

    //! Draws the leaves to file, without duplicating any lines
    void drawPrettyMesh(std::string fname);

    //! \copydoc drawPrettyMesh()
    void drawPrettyMesh(std::ostream& os);

#endif


//! Draws the leaves of the tree whose rank is the current node
void drawPartialMesh(std::string fname);

//! Draws the leaves of the tree within a range
void drawPartialMesh(std::string fname, cell_t from, cell_t to);
//! \copydoc drawPartialMesh(std::string, cell_t, cell_t)
void drawPartialMesh(std::ostream& os, cell_t from, cell_t to);


//! Draws the ghost cells.
void drawGhosts(std::string fname, int n);
//! \copydoc drawGhosts()
void drawGhosts(std::ostream& os, int n);


//! Draws the border cells
void drawBorder(std::string fname, int n);
//! \copydoc drawBorder()
void drawBorder(std::ostream& os, int n);


//! Draws the space filling curve of the tree
void drawCurve(std::string fname);
//! \copydoc drawCurve()
void drawCurve(std::ostream& os);

//! Draws the curve of the leaves in the tree whose rank is the current node
void drawPartialCurve(std::string fname);
//! \copydoc drawPartialCurve(std::string)
void drawPartialCurve(std::ostream& os);
//! Draws the curve of the leaves in the tree within a range
void drawPartialCurve(std::string fname, cell_t from, cell_t to);
//! \copydoc drawPartialCurve(std::string, cell_t, cell_t)
void drawPartialCurve(std::ostream& os, cell_t from, cell_t to);


//! Draws the boundary cells
void drawBoundaries(std::string fname);
//! Draws the boundary cells of a single boundary
void drawBoundary(std::ostream& os, int b);

//! Draws the normalised data to an image file
void drawMatrix(std::string fname, int imgW, int imgH, fnDataNorm dataNorm);
//! \copydoc drawMatrix
void drawMatrix(std::ostream& os, int imgW, int imgH, fnDataNorm dataNorm);

//! Plots the data
void plotMatrix(std::string fname, fnData dt);
//! \copydoc plotMatrix
void plotMatrix(std::ostream& os, fnData dt);


//! Draws the cells required to perform the poisson calculations
void drawPoissonNeighbourhood(std::string fname, cell_t cl);
//! \copydoc drawPoissonNeighbourhood
void drawPoissonNeighbourhood(std::ostream& os, cell_t cl);


} // namespace tftt

#endif
