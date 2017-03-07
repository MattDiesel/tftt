

#ifndef TFTT_FACEREF_H
#define TFTT_FACEREF_H


#include "tftt.h"
#include "treeface.h"


namespace tftt  {


struct FaceRef {
    TreeFace* face;

    FaceRef(TreeFace*);

    facedata_t* operator->();
    facedata_t const* operator->() const;

};


} // namespace tftt


#endif
