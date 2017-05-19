
#ifdef TFTT_FACES

#include "config.h"
#include "structure/treeface.h"

#include "faceref.h"


namespace tftt {

FaceRef::FaceRef(TreeFace* tf) : face(tf)
{
}


facedata_t* FaceRef::operator->()
{
    return &face->data;
}

facedata_t const* FaceRef::operator->() const
{
    return &face->data;
}


} // namespace tftt


#endif
