
#ifndef TFTT_ITER_LEAVES_H
#define TFTT_ITER_LEAVES_H


#include "../config.h"
#include "../cellref.h"
#include "cellrefiterator.h"


namespace tftt {


struct tagLeaves {
    class leaf_iterator : public cellref_iterator {
        void next();
    public:
        leaf_iterator(CellRef c);

        leaf_iterator operator++();
        leaf_iterator operator++(int junk);
    };

    leaf_iterator begin();
    leaf_iterator end();
};

extern tagLeaves leaves;


} // namespace tftt


#endif
