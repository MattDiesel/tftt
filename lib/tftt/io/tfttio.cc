
#include <string>
#include <vector>
#include <fstream>
#include <ostream>
#include <stdexcept>

#include <iostream> // Todo: Remove

#include "util/formatstring.h"

#include "../config.h"
#include "../structure/tree.h"
#include "../structure/treecell.h"
#include "../structure/treegroup.h"
#include "../fttcore.h"
#include "../iter/all.h"
#include "../gray.h"
#include "../tfttops.h"

#include "tfttio.h"


namespace tftt {


void drawMeshSub2d(std::ostream& os, TreeGroup* gr, double w, double h, double x, double y)
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
        drawMeshSub2d(os, gr->cells[0].children, w*0.5, h*0.5, x, y);
        drawMeshSub2d(os, gr->cells[1].children, w*0.5, h*0.5, x+w*0.5, y);
        drawMeshSub2d(os, gr->cells[2].children, w*0.5, h*0.5, x, y+h*0.5);
        drawMeshSub2d(os, gr->cells[3].children, w*0.5, h*0.5, x+w*0.5, y+h*0.5);
    }
}


void drawMesh(std::string fname)
{
    std::ofstream ofs(fname);
    drawMesh(ofs);
}

void drawMesh(std::ostream& os)
{
    if (DIM == 2) {
        drawMeshSub2d(os, gtree.root, gtree.size[0], gtree.size[1], 0.0, 0.0);
    }
    else {
        for (auto& c : leaves) {
            drawCell(os, c);
        }
    }
}


#if DIM == 2

void drawPrettyMesh(std::string fname)
{
    std::ofstream ofs(fname);
    drawPrettyMesh(ofs);
}

void drawPrettyMesh(std::ostream& os, cell_t cl)
{
    if (!cl.hasChildren()) return;

    os << cl.centre(0) << " " << cl.origin(1) << "\n"
       << cl.centre(0) << " " << (cl.origin(1)+cl.size(1)) << "\n\n"
       << cl.origin(0) << " " << cl.centre(1) << "\n"
       << (cl.origin(0)+cl.size(0)) << " " << cl.centre(1) << "\n\n";

    for (auto clch : *cl.children()) {
        drawPrettyMesh(os, clch);
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

void drawPrettyMesh(std::ostream& os)
{
    drawPrettyMesh(os, CellRef(-1));
}

#endif


void drawCell(std::ostream& os, cell_t const& c)
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


void drawPartialMesh(std::string fname, cell_t from, cell_t to)
{
    std::ofstream ofs(fname);
    drawPartialMesh(ofs, from, to);
}

void drawPartialMesh(std::string fname)
{
    std::ofstream ofs(fname);
    drawPartialMesh(ofs, gtree.firstActive, gtree.lastActive);
}

void drawPartialMesh(std::ostream& os, cell_t from, cell_t to)
{
    auto bgn = tagCurve::curve_iterator(from);
    auto end = tagCurve::curve_iterator(to);
    end++;

    for (auto& cl = bgn; cl != end; cl++) {
        drawCell(os, *cl);
    }
}


void drawGhosts(std::string fname, int n)
{
    std::ofstream ofs(fname);
    drawGhosts(ofs, n);
}

void drawGhosts(std::ostream& os, int n)
{
    if (!gtree.ghosts[n].empty()) {
        for (auto& c : gtree.ghosts[n]) {
            drawCell(os, c);
        }
    }
}


void drawBorder(std::string fname, int n)
{
    std::ofstream ofs(fname);
    drawBorder(ofs, n);
}

void drawBorder(std::ostream& os, int n)
{
    if (!gtree.borders[n].empty()) {
        for (auto& c : gtree.borders[n]) {
            drawCell(os, c);
        }
    }
}


void drawCurve(std::string fname)
{
    std::ofstream ofs(fname);
    drawCurve(ofs);
}

void drawCurve(std::ostream& os)
{
    for (auto& c : curve) {
        os << c.centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << c.centre(d);
        }

        os << "\n";
    }
}


void drawPartialCurve(std::string fname)
{
    std::ofstream ofs(fname);
    drawPartialCurve(ofs);
}

void drawPartialCurve(std::ostream& os)
{
    for (auto& c : activecurve) {
        os << c.centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << c.centre(d);
        }

        os << "\n";
    }
}

void drawPartialCurve(std::string fname, cell_t from, cell_t to)
{
    std::ofstream ofs(fname);
    drawPartialCurve(ofs, from, to);
}

void drawPartialCurve(std::ostream& os, cell_t from, cell_t to)
{
    auto bgn = tagCurve::curve_iterator(from);
    auto end = tagCurve::curve_iterator(to);
    end++;

    for (auto& cl = bgn; cl != end; cl++) {
        os << cl->centre(0);
        for (int d = 1; d < DIM; d++) {
            os << " " << cl->centre(d);
        }

        os << "\n";
    }
}


void drawBoundaries(std::string fname)
{
    std::ofstream ofs(fname);
    for (int b = 0; b < 2*DIM; b++) {
        drawBoundary(ofs, b);
    }
}

void drawBoundary(std::ostream& os, int b)
{
    for (auto& c : boundaryCells(b)) {
        drawCell(os, c);
    }
}


void drawMatrix(std::string fname, int imgW, int imgH, fnDataNorm dataNorm)
{
    std::ofstream ofs(fname, std::ios::binary);
    drawMatrix(ofs, imgW, imgH, dataNorm);
}

void drawMatrix(std::ostream& os, int imgW, int imgH, fnDataNorm dataNorm)
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
                bmp[y*imgW+x] = (unsigned char)(dataNorm(cell.data(), 255));
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


void plotMatrix(std::string fname, fnData dt)
{
    std::ofstream ofs(fname);
    plotMatrix(ofs, dt);
}

void plotMatrix(std::ostream& os, fnData dt)
{
    for (auto& cl : activecurve) {
        os << cl.centre(0) << " " << cl.centre(1) << " " << dt(cl.data()) << "\n";
    }
}


void drawPoissonNeighbourhood(std::string fname, cell_t cl)
{
    std::ofstream ofs(fname);
    drawPoissonNeighbourhood(ofs, cl);
}

void drawPoissonNeighbourhood(std::ostream& os, cell_t cl)
{

    TreeCell* tc = &cl.group->cells[cl.index];

    for (int n = 0; n < tc->poisNgbC; n++) {
        drawCell(os, tc->poisNgb[n]);
    }

}





} // namespace tftt
