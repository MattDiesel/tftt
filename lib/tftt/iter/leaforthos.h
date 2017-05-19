
#ifndef TFTT_ITER_LEAFORTHOS_H
#define TFTT_ITER_LEAFORTHOS_H


#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct tagLeafOrthos {
    class ortho_iterator {
        CellRef cr;
        CellRef ortho;

        void next();
    public:
        ortho_iterator(CellRef c);

        ortho_iterator operator++();
        ortho_iterator operator++(int junk);
        CellRef& operator*();
        CellRef* operator->();
        bool operator==(const ortho_iterator& rhs);
        bool operator!=(const ortho_iterator& rhs);
    };

    ortho_iterator begin();
    ortho_iterator end();
};

extern tagLeafOrthos leaforthos;


} // namespace tftt


#endif
