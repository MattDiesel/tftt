
#ifndef TFTT_VERTEXREF_H
#define TFTT_VERTEXREF_H

#ifdef TFTT_VERTICES


#include "config.h"
#include "structure/treevertex.h"


namespace tftt {


struct VertexRef {

    TreeVertex* vert;

    VertexRef(TreeVertex* tv);

};


} // namespace tftt


#endif
#endif
