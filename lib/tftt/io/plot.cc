
#include <string>
#include <fstream>
#include <ostream>

#include "../config.h"
#include "../cellref.h"
#include "../structure/treegroup.h"
#include "../structure/tree.h"
#include "../gray.h"
#include "../iter/all.h"

#include "plot.h"


namespace tftt {
namespace plot {


template<typename T>
void curveIter(std::ostream& os, T cells)
{
    for (auto&& c : cells) {
        os << c.centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << c.centre(d);
        }

        os << "\n";
    }
}


template<typename T>
void curveIter(std::string fname, T cells)
{
    std::ofstream ofs(fname);
    curveIter(ofs, cells);
}


template<typename IT = tagCurve::curve_iterator>
void curveRangeIter(std::ostream& os, cell_t from, cell_t to)
{
    auto bgn = IT(from);
    auto end = IT(to);

    if (to.isValid())
        end++;

    for (auto&& cl = bgn; cl != end; cl++) {
        os << cl->centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << cl->centre(d);
        }

        os << "\n";
    }
}


template<typename IT = tagCurve::curve_iterator>
void curveRangeIter(std::string fname, cell_t from, cell_t to)
{
    std::ofstream ofs(fname);
    curveRangeIter<IT>(ofs, from, to);
}


template<typename T>
void meshIter(std::ostream& os, T cells)
{
    for (auto& c : cells) {
        cellRect(os, c);
    }
}


template<typename T>
void meshIter(std::string fname, T cells)
{
    std::ofstream ofs(fname);
    meshIter(ofs, cells);
}


template<typename IT = tagCurve::curve_iterator>
void meshRangeIter(std::ostream& os, cell_t from, cell_t to)
{
    auto bgn = IT(from);
    auto end = IT(to);

    if (to.isValid())
        end++;

    for (auto& cl = bgn; cl != end; cl++) {
        cellRect(os, *cl);
    }
}


template<typename IT = tagCurve::curve_iterator>
void meshRangeIter(std::string fname, cell_t from, cell_t to)
{
    std::ofstream ofs(fname);
    curveRangeIter<IT>(ofs, from, to);
}


void cellRect(std::ostream& os, cell_t const& c)
{
    double vtc[1<<DIM][DIM];

    c.vertices(vtc);

    int vg;
    for (int v = 0; v < 1<<DIM; v++) {
        vg = utils::toGray(v);

        // os << c.vertex(utils::toGray(v), 0);
        os << vtc[vg][0];
        for (int d = 1; d < DIM; d++) {
            // os << " " << c.vertex(vg, d);
            os << ' ' << vtc[vg][d];
        }
        os << '\n';
    }

    // os << c.vertex(0, 0);
    os << vtc[0][0];
    for (int d = 1; d < DIM; d++) {
        // os << " " << c.vertex(0, d);
        os << ' ' << vtc[0][d];
    }
    os << "\n\n";
}


void mesh2d(std::ostream& os, TreeGroup* gr, double w, double h, double x, double y)
{
    // if (!gr) {
    os << x << ' ' << y << '\n'
       << (x+w) << ' ' << y << '\n'
       << (x+w) << ' ' << (y+h) << '\n'
       << x << ' ' << (y+h) << '\n'
       << x << ' ' << y << "\n\n";
    // }
    // else {
    if (gr) {
        mesh2d(os, gr->cells[0].children, w*0.5, h*0.5, x, y);
        mesh2d(os, gr->cells[1].children, w*0.5, h*0.5, x+w*0.5, y);
        mesh2d(os, gr->cells[2].children, w*0.5, h*0.5, x, y+h*0.5);
        mesh2d(os, gr->cells[3].children, w*0.5, h*0.5, x+w*0.5, y+h*0.5);
    }
}


void mesh(std::string fname)
{
    std::ofstream ofs(fname);
    mesh(ofs);
}


void mesh(std::ostream& os)
{
    if (DIM == 2) {
        mesh2d(os, gtree.root, gtree.size[0], gtree.size[1], 0.0, 0.0);
    }
    else {
        meshIter(os, leaves);
    }
}


#if DIM == 2

void prettyMesh(std::string fname)
{
    std::ofstream ofs(fname);
    prettyMesh(ofs);
}


void prettyMesh(std::ostream& os, cell_t cl)
{
    if (!cl.hasChildren()) return;

    os << cl.centre(0) << " " << cl.origin(1) << "\n"
       << cl.centre(0) << " " << (cl.origin(1)+cl.size(1)) << "\n\n"
       << cl.origin(0) << " " << cl.centre(1) << "\n"
       << (cl.origin(0)+cl.size(0)) << " " << cl.centre(1) << "\n\n";

    for (auto clch : *cl.children()) {
        prettyMesh(os, clch);
    }

    constexpr int sides[4][2] = {
        {0, 2},
        {1, 3},
        {0, 1},
        {2, 3}
    };

    if (cl.isRoot()) {
        // Draw borders
    }
    else {
        cell_t nbcl;
        for (int nb = 0; nb < (1<<DIM); nb++) {
            nbcl = cl.neighbour(nb);
            if (nbcl.parent() == cl.parent()) continue;

            if (nbcl.level() == cl.level() && !nbcl.hasChildren()) {
                os << cl.vertexPoint(sides[nb][0], 0) << " " << cl.vertexPoint(sides[nb][0], 1) << "\n"
                   << cl.vertexPoint(sides[nb][1], 0) << " " << cl.vertexPoint(sides[nb][1], 1) << "\n\n";
            }
        }
    }
}


void prettyMesh(std::ostream& os)
{
    prettyMesh(os, CellRef(-1));
}

#endif


void hilbert(std::string fname)
{
    std::ofstream ofs(fname);
    hilbert(ofs);
}


void hilbert(std::ostream& os)
{
    curveIter(os, curve);
}


void boundariesMesh(std::string fname)
{
    std::ofstream ofs(fname);
    for (int b = 0; b < 2*DIM; b++) {
        boundaryMesh(ofs, b);
    }
}

void boundaryMesh(std::ostream& os, int b)
{
    meshIter(os, boundaryCells(b));
}


void poissonNeighboursMesh(std::string fname, cell_t cl)
{
    std::ofstream ofs(fname);
    poissonNeighboursMesh(ofs, cl);
}


void poissonNeighboursMesh(std::ostream& os, cell_t cl)
{
    meshIter(os, cl.poissonNeighbourhood());
}


void partialMesh(std::string fname, cell_t from, cell_t to)
{
    std::ofstream ofs(fname);
    partialMesh(ofs, from, to);
}


void partialMesh(std::string fname)
{
    std::ofstream ofs(fname);
    meshIter(ofs, activecurve);
}


void partialMesh(std::ostream& os, cell_t from, cell_t to)
{
    meshRangeIter(os, from, to);
}


void partialHilbert(std::string fname)
{
    std::ofstream ofs(fname);
    partialHilbert(ofs);
}


void partialHilbert(std::ostream& os)
{
    curveIter(os, activecurve);
}


void partialHilbert(std::string fname, cell_t from, cell_t to)
{
    std::ofstream ofs(fname);
    partialHilbert(ofs, from, to);
}


void partialHilbert(std::ostream& os, cell_t from, cell_t to)
{
    curveRangeIter(os, from, to);
}


void ghostMesh(std::string fname, int n)
{
    std::ofstream ofs(fname);
    ghostMesh(ofs, n);
}


void ghostMesh(std::ostream& os, int n)
{
    if (!gtree.ghosts[n].empty())
        meshIter(os, gtree.ghosts[n]);
}


void borderMesh(std::string fname, int n)
{
    std::ofstream ofs(fname);
    borderMesh(ofs, n);
}


void borderMesh(std::ostream& os, int n)
{
    if (!gtree.borders[n].empty())
        meshIter(os, gtree.borders[n]);
}


} // namespace plot
} // namespace tftt
