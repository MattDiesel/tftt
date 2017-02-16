
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


struct tagBoundaryLeaves {
    int b;

    class bleaf_iterator {
        CellRef cr;
        int b;


        void next();
    public:
        bleaf_iterator(CellRef c, int bnd);

        bleaf_iterator operator++();
        bleaf_iterator operator++(int junk);
        CellRef& operator*();
        CellRef* operator->();
        bool operator==(const bleaf_iterator& rhs);
        bool operator!=(const bleaf_iterator& rhs);
    };

    bleaf_iterator begin();
    bleaf_iterator end();
};

tagBoundaryLeaves boundaryCells(int b);





struct tagCurve {
    class curve_iterator {
        CellRef cr;

        void next();
    public:
        curve_iterator(CellRef c);

        curve_iterator operator++();
        curve_iterator operator++(int junk);
        CellRef& operator*();
        CellRef* operator->();
        bool operator==(const curve_iterator& rhs);
        bool operator!=(const curve_iterator& rhs);
    };

    curve_iterator begin();
    curve_iterator end();
};

extern tagCurve curve;




} // namespace tftt


#endif
