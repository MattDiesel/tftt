
#ifndef TFTT_ITER_NEIGHBOURS_H
#define TFTT_ITER_NEIGHBOURS_H


#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct tagNeighbours {
    cell_t cr;

    tagNeighbours(cell_t c);

    class neighb_iterator {
        cell_t cl;
        int nb;
        cell_t cr;

        void next();

    public:
        neighb_iterator(cell_t c, int n);

        neighb_iterator operator++();
        neighb_iterator operator++(int junk);
        cell_t& operator*();
        cell_t* operator->();
        bool operator==(const neighb_iterator& rhs);
        bool operator!=(const neighb_iterator& rhs);
    };

    neighb_iterator begin();
    neighb_iterator end();
};


} // namespace tftt


#endif
