

#ifndef TFTT_FACES_H
#define TFTT_FACES_H


#ifndef TFTT_LEAVES
#define TFTT_LEAVES


#include "tftt.h"


namespace tftt {


struct tagFaces {
    class face_iterator {
        CellRef cr;

        void next();
    public:
        face_iterator(CellRef c);

        face_iterator operator++();
        face_iterator operator++(int junk);
        CellRef& operator*();
        CellRef* operator->();
        bool operator==(const face_iterator& rhs);
        bool operator!=(const face_iterator& rhs);
    };

    face_iterator begin();
    face_iterator end();
};

extern tagFaces faces;


} // namespace tftt


#endif


#endif
