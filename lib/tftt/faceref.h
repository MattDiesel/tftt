

#ifndef TFTT_FACEREF_H
#define TFTT_FACEREF_H

#ifdef TFTT_FACES

#include "config.h"
#include "structure/treeface.h"


namespace tftt  {


struct FaceRef {
    TreeFace* face;

    FaceRef(TreeFace*);

    facedata_t* operator->();
    facedata_t const* operator->() const;

};


} // namespace tftt


#endif
#endif
