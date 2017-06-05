
#include <string>
#include <ostream>
#include <fstream>

#include "util/formatstring.h"

#include "../config.h"
#include "../cellref.h"
#include "../structure/tree.h"
#include "../iter/all.h"

#include "tikz.h"


namespace tftt {
namespace tikz {


void meshLayer(std::ostream& os, cell_t cl, int layer)
{
    if (cl.level() == layer) {
        os << "\\node at ("
           << cl.centre(0) << "," << cl.centre(1) << ") {};\n";
    }
    else if (cl.hasChildren()) {
        for (auto ch : *cl.children()) {
            meshLayer(os, ch, layer);
        }
    }
}


void meshLayer(std::ostream& ofs, int layer)
{
    cell_t root = cell_t(-1);

    double w = gtree.size[0] / (2 << layer);
    double h = gtree.size[1] / (2 << layer);

    // ofs << "\\begin{scope}[every node/.append style={minimum width=" << w
    //     << "cm-\\the\\pgflinewidth,minimum height=" << h
    //     << "cm-\\the\\pgflinewidth,inner sep=0,draw,shape=rectangle,anchor=center}]\n";

    ofs << "\\begin{scope}[every node/.append style={minimum width=" << w
        << "cm,minimum height=" << h
        << "cm,inner sep=0,draw,shape=rectangle,anchor=center}]\n";

    meshLayer(ofs, root, layer);

    ofs << "\\end{scope}\n";
}


void meshLayer(std::string fname, int layer)
{
    std::ofstream ofs(fname);
    meshLayer(ofs, layer);
}


void mesh(std::string fname, int maxDepth)
{
    std::ofstream ofs(fname);

    for (int d = 0; d < maxDepth; d++) {
        meshLayer(ofs, d);
    }
}


template<class T>
void curve(std::string fname, T set)
{
    std::ofstream ofs(fname);

    cell_t prev;
    for (auto cl : set) {
        if (!prev.isValid()) {
            prev = cl;
            continue;
        }

        ofs << "\\draw (" << prev.centre(0) << "," << prev.centre(1)
            << ") -- (" << cl.centre(0) << "," << cl.centre(1) << ");\n";

        prev = cl;
    }

    ofs.close();
}


void hilbert(std::string fname)
{
    curve(fname, tftt::curve);
}


void morton(std::string fname)
{
    curve(fname, tftt::leaves);
}




void meshLayer_if(std::ostream& os, cell_t cl, int layer, fnCheckCell pred)
{
    if (cl.level() == layer) {
        if (pred(cl)) {
            os << "\\node at ("
               << cl.centre(0) << "," << cl.centre(1) << ") {};\n";
        }
    }
    else if (cl.hasChildren()) {
        for (auto ch : *cl.children()) {
            meshLayer_if(os, ch, layer, pred);
        }
    }
}


void meshLayer_if(std::ostream& ofs, int layer, fnCheckCell pred)
{
    cell_t root = cell_t(-1);

    double w = gtree.size[0] / (2 << layer);
    double h = gtree.size[1] / (2 << layer);

    // ofs << "\\begin{scope}[every node/.append style={minimum width=" << w
    //     << "cm-\\the\\pgflinewidth,minimum height=" << h
    //     << "cm-\\the\\pgflinewidth,inner sep=0,draw,shape=rectangle,anchor=center}]\n";

    ofs << "\\begin{scope}[every node/.append style={minimum width=" << w
        << "cm,minimum height=" << h
        << "cm,inner sep=0,draw,shape=rectangle,anchor=center}]\n";

    meshLayer_if(ofs, root, layer, pred);

    ofs << "\\end{scope}\n";
}


void meshLayer_if(std::string fname, int layer, fnCheckCell pred)
{
    std::ofstream ofs(fname);
    meshLayer_if(ofs, layer, pred);
}


void mesh_if(std::string fname, int maxDepth, fnCheckCell pred)
{
    std::ofstream ofs(fname);

    for (int d = 0; d < maxDepth; d++) {
        meshLayer_if(ofs, d, pred);
    }
}


} // namespace tikz
} // namespace tftt
