
#ifndef TFTT_ITER_CELLREF_H
#define TFTT_ITER_CELLREF_H


#include "../config.h"
#include "../cellref.h"


namespace tftt {


class cellref_iterator {
protected:
    cell_t cr;
public:
    cellref_iterator(cell_t c);

    cell_t& operator*();
    cell_t* operator->();
    bool operator==(const cellref_iterator& rhs);
    bool operator!=(const cellref_iterator& rhs);
};



} // namespace tftt


#endif
