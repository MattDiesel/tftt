
#ifndef TFTT_ITER_CELLREF_H
#define TFTT_ITER_CELLREF_H


#include "../config.h"
#include "../cellref.h"


namespace tftt {


class cellref_iterator {
protected:
    CellRef cr;
public:
    cellref_iterator(CellRef c);

    CellRef& operator*();
    CellRef* operator->();
    bool operator==(const cellref_iterator& rhs);
    bool operator!=(const cellref_iterator& rhs);
};



} // namespace tftt


#endif
