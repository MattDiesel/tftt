
#ifndef TFTT_ITER_NEIGHBOURS_H
#define TFTT_ITER_NEIGHBOURS_H


#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct tagNeighbours {
    CellRef cr;

    tagNeighbours(CellRef c);

    class neighb_iterator {
        CellRef cl;
        int nb;
        CellRef cr;

        void next();

    public:
        neighb_iterator(CellRef c, int n);

        neighb_iterator operator++();
        neighb_iterator operator++(int junk);
        CellRef& operator*();
        CellRef* operator->();
        bool operator==(const neighb_iterator& rhs);
        bool operator!=(const neighb_iterator& rhs);
    };

    neighb_iterator begin();
    neighb_iterator end();
};


} // namespace tftt


#endif
