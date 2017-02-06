
#ifndef TFTT_LEAVES
#define TFTT_LEAVES


#include "tftt.h"


namespace tftt {


struct tagLeaves {
    class leaf_iterator {
        CellRef cr;

        void next();
    public:
        leaf_iterator(CellRef c);

        leaf_iterator operator++();
        leaf_iterator operator++(int junk);
        CellRef& operator*();
        CellRef* operator->();
        bool operator==(const leaf_iterator& rhs);
        bool operator!=(const leaf_iterator& rhs);
    };

    leaf_iterator begin();
    leaf_iterator end();
};

extern tagLeaves leaves;


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
