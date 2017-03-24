
#include <iostream>
#include <sstream>

#define DIM 2
#include "tftt/treeid.h"

typedef tftt::TreeId<uint64_t> ident_t;


template<typename T>
bool tfetch(std::string val, T& ret)
{
    std::istringstream in(val);
    in >> ret;
    return !in.fail();
}


int main(int argc, char const* argv[])
{

    if (argc < 3) {
        std::cout << "Usage: pos2id level x y\n";
        return 1;
    }

    int level;
    double x, y;

    tfetch(argv[1], level);
    tfetch(argv[2], x);
    tfetch(argv[3], y);

    std::cout << "Examining cell at {" << x << "," << y << "} depth " << level << "\n";

    double xl = 0.0, xu = 1.0, yl = 0.0, yu = 1.0;
    int depth = 0;
    ident_t id;
    int ch;

    while (depth <= level) {
        ch = 0;

        if (x < (xl+xu)*0.5)
            xu = (xl+xu)*0.5;
        else {
            xl = (xl+xu)*0.5;
            ch += 1;
        }

        if (y < (yl+yu)*0.5)
            yu = (yl+yu)*0.5;
        else {
            yl = (yl+yu)*0.5;
            ch += 2;
        }

        if (depth)
            id = id.child(ch);
        else
            id.id |= ch;

        depth++;
    }

    std::cout << "Cell " << id.id << '\n';
    std::cout << "\tOrigin: " << xl << ", " << yl << '\n';
    std::cout << "\tCentre: " << (xl+xu)*0.5 << ", " << (yl+yu)*0.5 << '\n';
    std::cout << "\tSize: " << (xu-xl) << " x " << (yu - yl) << '\n';

    return 0;
}

