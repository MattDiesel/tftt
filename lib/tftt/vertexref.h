
#ifndef TFTT_VERTEXREF_H
#define TFTT_VERTEXREF_H


#include "tftt.h"
#include "treevertex.h"


namespace tftt {


struct VertexRef {

    TreeVertex* vert;

    VertexRef(TreeVertex* tv);

};


} // namespace tftt


#endif
