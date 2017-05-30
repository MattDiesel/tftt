
#ifndef TFTT_ITER_CURVE_H
#define TFTT_ITER_CURVE_H


#include "../config.h"
#include "../cellref.h"
#include "cellrefiterator.h"


namespace tftt {


struct tagCurve {
    class curve_iterator : public cellref_iterator {
        void next();
        void prev();
    public:
        curve_iterator(cell_t c);

        curve_iterator operator++();
        curve_iterator operator++(int junk);

        curve_iterator operator--();
        curve_iterator operator--(int junk);
    };

    curve_iterator begin();
    curve_iterator end();

    curve_iterator rbegin();
    curve_iterator rend();
};

extern tagCurve curve;


} // namespace tftt


#endif
