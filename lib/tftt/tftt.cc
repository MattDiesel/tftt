
#include <fstream>
#include <stdexcept>
#include <set>
#include <cmath>
#include <iostream> // TODO: Remove

#include "util/formatstring.h"

#include "tftt.h"
#include "structure/tree.h"
#include "structure/treegroup.h"
#include "cellref.h"
#include "gray.h"
#include "hilbert.h"


namespace tftt {


bool checkAround(cell_t cl, int dist, fnCheckCell check)
{
    cell_t tmp;
    for (int nb = 0; nb < 2*DIM; nb++) {
        tmp = cl;
        for (int p = 0; p < dist; p++) {
            tmp = tmp.neighbour(nb);
            if (tmp.isBoundary()) break;
            if (!check(tmp)) return false;
        }
    }

    return true;
}


bool findAround(cell_t cl, int dist, fnCheckCell check)
{
    cell_t tmp;
    for (int nb = 0; nb < 2*DIM; nb++) {
        tmp = cl;
        for (int p = 0; p < dist; p++) {
            tmp = tmp.neighbour(nb);
            if (tmp.isBoundary()) break;
            if (check(tmp)) return true;
        }
    }

    return false;
}


double interpFace(cell_t cl, int fc, fnData dt)
{
    cell_t ngb = cl.neighbour(fc);
    int clLvl = cl.level();
    int ngbLvl = ngb.level();
    double ret;

    if (ngb.hasChildren()) {
        // Neighbour more refined
        // Average children on that face
        // Todo: Generalise to 3d
        // std::cout << "2\n";

        ret = 0.0;
        for (int n = 0; n < 2; n++) {
            ret += dt(ngb.childOnFace(fc ^ 1, n).data()) / 3.0;
        }

        ret += dt(cl.data()) / 3.0;
    }
    else if (clLvl == ngbLvl) {
        // Same level
        // std::cout << "1\n";

        ret = 0.5*(dt(cl.data()) + dt(ngb.data()));
    }
    else {
        // Neighbour less refined

        int awayDir = (fc ^ 2) & ~1;
        if (cl.index() & (1<<(awayDir >> 1))) awayDir ^= 1;

        cell_t ngbngb = ngb.neighbour(awayDir);

        ret = 0.0;
        if (ngbngb.hasChildren()) {
            // std::cout << "4\n";

            for (int n = 0; n < 2; n++) {
                ret += dt(ngbngb.childOnFace(awayDir ^ 1, n).data()) / 18.0;
            }
            ret += dt(ngb.data()) * 2.0 / 9.0;
            ret += dt(cl.data()) * 2.0 / 3.0;
        }
        else {
            // std::cout << "3\n";

            ret += dt(ngbngb.data()) / 12.0;
            ret += dt(ngb.data()) / 4.0;
            ret += dt(cl.data()) * 2.0 / 3.0;
        }
    }

    return ret;
}


double interpChild(cell_t cl, int ch, int forDir, fnData dt)
{
    cell_t forNbCl = cl.neighbour(forDir);

    // if (!forNbCl.hasChildren()) {
    //     return cl.data();
    // }

    forNbCl = forNbCl.child(ch ^ (1 << (forDir >> 1)));

    int awayDir;
    cell_t awayNb;

    double ret = dt(cl.data());
    double intDir;

    for (int d = 0; d < DIM; d++) {
        if (d == forDir >> 1) continue; // Interpolate in desired direction last

        awayDir = (d << 1);
        if (ch & awayDir) awayDir++;

        awayNb = cl.neighbour(awayDir);

        if (awayNb.isBoundary()) {
            continue; // Will just interpolate to the same.
        }

        if (awayNb.hasChildren()) {
            // Interpolate between boundary children
            // Todo: Generalise
            int c1 = ch ^ (1 << d);
            int c2 = c1 ^ (1 << (forDir >> 1));
            intDir = (dt(awayNb.child(c1).data())
                      + dt(awayNb.child(c2).data()))*0.5;

            ret += (intDir - ret) / 3.0;
        }
        else {
            intDir = dt(awayNb.data());
            ret += (intDir - ret) / 4.0;
        }
    }

    // Finally interpolate in desired direction
    intDir = dt(forNbCl.data());
    ret += (intDir - ret) / 3.0;

    return ret;
}


double interpALEVertex(cell_t cl, int v, fnData dt)
{
    cell_t ngb = cl.neighbour(v | 1);
    cell_t c2;
    if (ngb.hasChildren()) {
        // Neighbour more refined, use child vertex
        ngb = ngb.child(0);
        return dt(ngb.data());
    }
    else if (ngb.level() < cl.level()) {
        // Neighbour less refined, average two vertices
        // Todo: Generalise to 3d
        c2 = ngb.neighbour(v ^ 2);
        if (c2.hasChildren())
            c2 = c2.child(0);

        return (dt(ngb.data()) + dt(c2.data()))*0.5;
    }
    else {
        // Same level.
        return dt(ngb.data());
    }
}


} // namespace tftt
