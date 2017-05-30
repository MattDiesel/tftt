
#ifndef TFTT_ITER_BOUNDARYCELLS_H
#define TFTT_ITER_BOUNDARYCELLS_H


#include "../config.h"
#include "../cellref.h"
#include "cellrefiterator.h"


namespace tftt {


struct tagBoundaryLeaves {
    int b;

    class bleaf_iterator : public cellref_iterator {
        int b;

        void next();
    public:
        bleaf_iterator(cell_t c, int bnd);

        bleaf_iterator operator++();
        bleaf_iterator operator++(int junk);
    };

    bleaf_iterator begin();
    bleaf_iterator end();
};

tagBoundaryLeaves boundaryCells(int b);


} // namespace tftt


#endif
