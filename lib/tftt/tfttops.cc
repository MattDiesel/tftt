

#include <cmath>
#include <stdexcept>

#include "config.h"
#include "treeid.h"
#include "cellref.h"
#include "structure/tree.h"
#include "structure/treegroup.h"
#include "structure/treecell.h"
#include "fttcore.h"
#include "adapt.h"
#include "iter/leaves.h"

#include "tfttops.h"


namespace tftt {


void calcFaceCoef(cell_t cl, TreeCell* tc, int fc);

void calcFaceCoefs(cell_t cl)
{
    TreeCell* tc = &cl.group->cells[cl.index];

    tc->poisNgbC = 0;

    for (int dir = 0; dir < 4; dir++) {
        calcFaceCoef(cl, tc, dir);
    }

    // Calc rhs
    double alphas = 0.0;
    for (int d = 0; d < 2*DIM; d++) {
        alphas += tc->poisAlpha[d]; // Todo: Rho?
        // std::cout << "\talpha[" << d << "] = " << tc.poisAlpha[d] << "\n";
    }

    tc->cenCoef = alphas;
}


void calcFaceCoef(cell_t cl, TreeCell* tc, int fc)
{
    cell_t ngb = cl.neighbour(fc);
    int clLvl = cl.level();
    int ngbLvl = ngb.level();

    cell_t tmp;

    double dbound = 1.0 / cl.size((fc >> 1) ^ 1);
    dbound *= dbound;

    if (ngb.isBoundary()) {
        if (gtree.isNeuman[fc]) {
            if (false) {
                // Fake Neuman
                tc->poisAlpha[fc] = 2.0*dbound;
                tc->poisNgb[tc->poisNgbC] = cl;
                tc->poisNgbDir[tc->poisNgbC] = fc & 1;
                tc->poisCoef[tc->poisNgbC++] = 2.0 * dbound;
            }
            else {
                tc->poisAlpha[fc] = 0.0;
            }
        }
        else {
            tc->poisAlpha[fc] = 2.0*dbound;
            tc->poisNgb[tc->poisNgbC] = ngb;
            tc->poisNgbDir[tc->poisNgbC] = fc & 1;
            tc->poisCoef[tc->poisNgbC++] = 2.0 * dbound;
        }

        return;
    }
    else if (ngb.hasChildren()) {
        // Neighbour more refined
        // Average children on that face
        // Todo: Generalise to 3d
        // std::cout << "2\n";

        if (false) {
            // Using summation of gradients
            tc->poisAlpha[fc] = 0.0;
            dbound *= 0.5;

            for (int d = 0; d < 2; d++) {
                cell_t ps = ngb.childOnFace(fc ^ 1, d);

                int awayDir = (fc ^ 2) & ~1;
                if (ps.index & (1<<(awayDir >> 1))) awayDir ^= 1;

                cell_t ngbngb = cl.neighbour(awayDir);

                if (ngbngb.hasChildren()) {
                    // std::cout << "2.4\n";

                    tc->poisAlpha[fc] += 8.0/9.0 * dbound;

                    for (int n = 0; n < 2; n++) {
                        tmp = ngbngb.childOnFace(awayDir ^ 1, n);

                        tc->poisNgb[tc->poisNgbC] = tmp;
                        tc->poisNgbDir[tc->poisNgbC] = fc & 1;
                        tc->poisCoef[tc->poisNgbC++] = (2.0/9.0) * dbound;

                        // std::cout << "Use " << tmp << " * " << ((1.0/18.0)*dbound) << "\n";
                    }

                    tc->poisNgb[tc->poisNgbC] = ps;
                    tc->poisNgbDir[tc->poisNgbC] = fc & 1;
                    tc->poisCoef[tc->poisNgbC++] = (4.0 / 3.0) * dbound;

                    // std::cout << "Use " << tmp << " * " << (2.0/9.0*dbound) << "\n";
                }
                else {
                    // std::cout << "2.3\n";

                    tc->poisAlpha[fc] += dbound;

                    tc->poisNgb[tc->poisNgbC] = ngbngb;
                    tc->poisNgbDir[tc->poisNgbC] = fc & 1;
                    tc->poisCoef[tc->poisNgbC++] = (1.0 / 3.0) * dbound;

                    // std::cout << "Use " << ngbngb << " * " << (1.0/12.0*dbound) << "\n";

                    tc->poisNgb[tc->poisNgbC] = ps;
                    tc->poisNgbDir[tc->poisNgbC] = fc & 1;
                    tc->poisCoef[tc->poisNgbC++] = (4.0 / 3.0) * dbound;

                    // std::cout << "Use " << ngbngb << " * " << (0.25*dbound) << "\n";

                }

            }
        }
        else {
            // Use interpolation of neighbours

            tc->poisAlpha[fc] = 4.0/3.0 * dbound;

            for (int n = 0; n < 2; n++) {
                tmp = ngb.childOnFace(fc ^ 1, n);

                tc->poisNgb[tc->poisNgbC] = tmp;
                tc->poisNgbDir[tc->poisNgbC] = fc & 1;
                tc->poisCoef[tc->poisNgbC++] = 2.0/3.0 * dbound;

                // std::cout << "Use " << tmp << " * " << (0.3333*dbound) << "\n";
            }
        }
    }
    else if (clLvl == ngbLvl) {
        // Same level
        // std::cout << "1\n";

        tc->poisAlpha[fc] = 1.0 * dbound;
        tc->poisNgb[tc->poisNgbC] = ngb;
        tc->poisNgbDir[tc->poisNgbC] = fc & 1;
        tc->poisCoef[tc->poisNgbC++] = 1.0 * dbound;

        // std::cout << "Use " << ngb << " * " << (0.5*dbound) << "\n";
    }
    else {
        // Neighbour less refined

        int awayDir = (fc ^ 2) & ~1;
        if (cl.index & (1<<(awayDir >> 1))) awayDir ^= 1;

        cell_t ngbngb = ngb.neighbour(awayDir);

        if (ngbngb.hasChildren()) {
            // std::cout << "4\n";

            tc->poisAlpha[fc] = 2.0/3.0 * dbound;

            for (int n = 0; n < 2; n++) {
                tmp = ngbngb.childOnFace(awayDir ^ 1, n);

                tc->poisNgb[tc->poisNgbC] = tmp;
                tc->poisNgbDir[tc->poisNgbC] = fc & 1;
                tc->poisCoef[tc->poisNgbC++] = (1.0/9.0) * dbound;

                // std::cout << "Use " << tmp << " * " << ((1.0/18.0)*dbound) << "\n";
            }

            tc->poisNgb[tc->poisNgbC] = ngb;
            tc->poisNgbDir[tc->poisNgbC] = fc & 1;
            tc->poisCoef[tc->poisNgbC++] = (4.0 / 9.0) * dbound;

            // std::cout << "Use " << tmp << " * " << (2.0/9.0*dbound) << "\n";
        }
        else {
            // std::cout << "3\n";

            tc->poisAlpha[fc] = 2.0/3.0 * dbound;

            tc->poisNgb[tc->poisNgbC] = ngbngb;
            tc->poisNgbDir[tc->poisNgbC] = fc & 1;
            tc->poisCoef[tc->poisNgbC++] = (1.0 / 6.0) * dbound;

            // std::cout << "Use " << ngbngb << " * " << (1.0/12.0*dbound) << "\n";

            tc->poisNgb[tc->poisNgbC] = ngb;
            tc->poisNgbDir[tc->poisNgbC] = fc & 1;
            tc->poisCoef[tc->poisNgbC++] = 0.5 * dbound;

            // std::cout << "Use " << ngbngb << " * " << (0.25*dbound) << "\n";

        }
    }
}


void relax(double omega, fnDataRef datafn, fnCell cellfn)
{
    double alphas, betas;
    double temp;

    for (auto& cl : tftt::activecurve) {
        TreeCell& tc = cl.group->cells[cl.index];

        betas = 0.0;
        for (int n = 0; n < tc.poisNgbC; n++) {
            betas += datafn(tc.poisNgb[n].data()) * tc.poisCoef[n]; // * (tc.poisNgbDir ? 1.0 : -1.0); // Todo: Rho?
        }
        alphas = tc.cenCoef;

        temp = (betas - cellfn(cl)) / alphas;

        datafn(cl.data()) = omega*temp + (1.0-omega)*datafn(cl.data());
    }
}


double resid(fnDataRef datafn, fnCell cellfn)
{
    TreeCell* tc;

    double max = 0.0;
    cell_t maxat;

    for (auto& cl : activecurve) {
        tc = &cl.group->cells[cl.index];

        cl->res = -cellfn(cl) - tc->cenCoef * datafn(cl.data());

        for (int n = 0; n < tc->poisNgbC; n++) {
            cl->res += datafn(tc->poisNgb[n].data()) * tc->poisCoef[n];
        }

        if (std::fabs(cl->res) > std::fabs(max)) {
            max = cl->res;
            maxat = cl;
        }
    }

    return max;
}


} // namespace tftt
