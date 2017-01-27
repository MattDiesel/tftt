
#include "tftt.h"
#include "tree.h"


namespace tftt {


Tree gtree;


Tree::~Tree() {
    delete root;
}


} // namespace tftt
