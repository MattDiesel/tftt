
#include "tftt.h"
#include "tree.h"


namespace tftt {


Tree gtree;


Tree::~Tree() {
    destroying = true;
    delete root;
}


} // namespace tftt
