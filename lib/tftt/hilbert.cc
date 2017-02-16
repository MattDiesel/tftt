
#include "tftt.h"

#include "hilbert.h"


namespace tftt {


int hilbOrient(int Dp, int ch) {
    static constexpr int orient[4] = {
        0x0021, 0x1310, 0x3202, 0x2133
    };

    return (orient[Dp] >> (ch*4)) & 15;
}


int hilbChild(int Dp, int hch) {
    static constexpr int child[4] = {
        0x1320, 0x2310, 0x1023, 0x2013
    };

    return (child[Dp] >> (hch*4)) & 15;
}


int hilbInvChild(int Dp, int ch) {
    static constexpr int child[4] = {
        0x2130, 0x2310, 0x0132, 0x0312
    };

    return (child[Dp] >> (ch*4)) & 15;
}


bool hilbIsLast(int Dp, int ch) {
    return hilbInvChild(Dp, ch) == (1<<DIM)-1;
}


bool hilbIsFirst(int Dp, int ch) {
    return hilbInvChild(Dp, ch) == 0;
}


} // namespace tftt
