
#ifndef TFTT_ITER_LEAFORTHOS_H
#define TFTT_ITER_LEAFORTHOS_H


#include "../config.h"
#include "../cellref.h"


namespace tftt {


struct tagLeafOrthos {
    class ortho_iterator {
        cell_t cr;
        cell_t ortho;

        void next();
    public:
        ortho_iterator(cell_t c);

        ortho_iterator operator++();
        ortho_iterator operator++(int junk);
        cell_t& operator*();
        cell_t* operator->();
        bool operator==(const ortho_iterator& rhs);
        bool operator!=(const ortho_iterator& rhs);
    };

    ortho_iterator begin();
    ortho_iterator end();
};

extern tagLeafOrthos leaforthos;


} // namespace tftt


#endif
