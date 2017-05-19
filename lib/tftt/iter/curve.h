
#ifndef TFTT_ITER_CURVE_H
#define TFTT_ITER_CURVE_H


#include "../config.h"
#include "../cellref.h"
#include "cellrefiterator.h"


namespace tftt {


struct tagCurve {
    class curve_iterator : public cellref_iterator {
        void next();
    public:
        curve_iterator(CellRef c);

        curve_iterator operator++();
        curve_iterator operator++(int junk);
    };

    curve_iterator begin();
    curve_iterator end();
};

extern tagCurve curve;


} // namespace tftt


#endif
