
#ifndef TFTT_IO_PLOT_H
#define TFTT_IO_PLOT_H


#include <string>
#include <ostream>

#include "../config.h"
#include "../cellref.h"


namespace tftt {
namespace plot {


//! Writes the coordinates of each cell vertex to a stream
//! \details When opened in gnuplot, this will outline the cell.
void cellRect(std::ostream& os, cell_t const& c);


//! Draws the leaves in the tree to a file
void mesh(std::string fname);
//! \copydoc mesh()
void mesh(std::ostream& os);


#if DIM == 2

    //! Draws the leaves to file, without duplicating any lines
    void prettyMesh(std::string fname);

    //! \copydoc prettyMesh()
    void prettyMesh(std::ostream& os);

#endif


//! Draws the space filling curve of the tree
void hilbert(std::string fname);
//! \copydoc hilbert()
void hilbert(std::ostream& os);


//! Draws the boundary cells
void boundariesMesh(std::string fname);
//! Draws the boundary cells of a single boundary
void boundaryMesh(std::ostream& os, int b);


//! Draws the cells required to perform the poisson calculations
void poissonNeighboursMesh(std::string fname, cell_t cl);
//! \copydoc poissonNeighboursMesh
void poissonNeighboursMesh(std::ostream& os, cell_t cl);


//! Draws the leaves of the tree whose rank is the current node
void partialMesh(std::string fname);
//! Draws the leaves of the tree within a range
void partialMesh(std::string fname, cell_t from, cell_t to);
//! \copydoc partialMesh(std::string, cell_t, cell_t)
void partialMesh(std::ostream& os, cell_t from, cell_t to);


//! Draws the curve of the leaves in the tree whose rank is the current node
void partialHilbert(std::string fname);
//! \copydoc partialHilbert(std::string)
void partialHilbert(std::ostream& os);
//! Draws the curve of the leaves in the tree within a range
void partialHilbert(std::string fname, cell_t from, cell_t to);
//! \copydoc partialHilbert(std::string, cell_t, cell_t)
void partialHilbert(std::ostream& os, cell_t from, cell_t to);


//! Draws the ghost cells.
void ghostMesh(std::string fname, int n);
//! \copydoc ghostMesh()
void ghostMesh(std::ostream& os, int n);


//! Draws the border cells
void borderMesh(std::string fname, int n);
//! \copydoc borderMesh()
void borderMesh(std::ostream& os, int n);


} // namespace plot
} // namespace tftt


#endif
