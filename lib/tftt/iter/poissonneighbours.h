
#ifndef TFTT_ITER_POISSONNEIGHBOURS_H
#define TFTT_ITER_POISSONNEIGHBOURS_H


#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct tagPoissonNeighbours {
    cell_t cr;

    tagPoissonNeighbours(cell_t c);

    class poisneighb_iterator {
        cell_t cl;
        int nb;
        cell_t cr;

        void next();

    public:
        poisneighb_iterator(cell_t c, int n);

        poisneighb_iterator operator++();
        poisneighb_iterator operator++(int junk);
        cell_t& operator*();
        cell_t* operator->();
        bool operator==(const poisneighb_iterator& rhs);
        bool operator!=(const poisneighb_iterator& rhs);
    };

    poisneighb_iterator begin();
    poisneighb_iterator end();
};


} // namespace tftt


#endif
