
#ifndef TFTT_ITER_POISSONNEIGHBOURS_H
#define TFTT_ITER_POISSONNEIGHBOURS_H


#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct tagPoissonNeighbours {
    CellRef cr;

    tagPoissonNeighbours(CellRef c);

    class poisneighb_iterator {
        CellRef cl;
        int nb;
        CellRef cr;

        void next();

    public:
        poisneighb_iterator(CellRef c, int n);

        poisneighb_iterator operator++();
        poisneighb_iterator operator++(int junk);
        CellRef& operator*();
        CellRef* operator->();
        bool operator==(const poisneighb_iterator& rhs);
        bool operator!=(const poisneighb_iterator& rhs);
    };

    poisneighb_iterator begin();
    poisneighb_iterator end();
};


} // namespace tftt


#endif
