

#include <cmath>
#include <stdexcept>
#include <iostream> // TODO: Remove

#include "tftt.h"
#include "tree.h"
#include "treegroup.h"

#include "tfttops.h"


namespace tftt {


cell_t find(ident_t idt)
{
    TreeGroup* gr = gtree.root;

    int n = 0;
    while (gr) {
        if (gr->id == idt.group()) {
            return CellRef(gr, (int)idt.orthant());
        }

        gr = gr->cells[idt.orthant(n++)].children;
    }

    return CellRef();
}

cell_t insert(ident_t idt)
{
    cell_t ret = CellRef(-1);

    while (ret.id().id != idt.id) {
        if (idt.level() <= ret.level())
            throw std::runtime_error("Badly formatted ID.");

        if (!ret.hasChildren()) {
            twoToOne(ret);
            refine(ret);
        }

        // std::cout << idt.orthant(ret.level() + 1) << std::endl;

        ret = ret.child(idt.orthant(ret.level() + 1));
    }

    return ret;
}


cell_t findmax(fnData dt, double* maxValRet)
{
    cell_t max;
    double maxVal = 0.0;
    double val;
    for (auto& cl : leaves) {
        val = dt(cl.data());
        if (val > maxVal) {
            max = cl;
            maxVal = val;
        }
    }

    if (maxValRet) *maxValRet = maxVal;
    return max;
}

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

    cl->cenCoef = alphas;
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
        if (options.isNeuman) {
            if (true) {
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


cell_t max(fnData dfn)
{
    double ret = 0.0, t;
    cell_t retc;
    for (auto& cl : tftt::leaves) {
        if (!retc.isValid()) {
            ret = dfn(cl.data());
            retc = cl;
        }
        else {
            t = dfn(cl.data());
            if (ret < t) {
                ret = t;
                retc = cl;
            }
        }
    }

    return retc;
}


void relax(double omega, fnDataRef datafn, fnCell cellfn)
{

    uint64_t target = 432345564227570477;

    double alphas, betas;
    double temp;

    for (auto& cl : tftt::leaves) {
        TreeCell& tc = cl.group->cells[cl.index];

        if (cl.id().id == target) {
            std::cout << "Relaxing target cell " << cl << "\n";
            std::cout << "\tdx = " << cl.size(0) << ", dy = " << cl.size(1) << "\n";
        }

        betas = 0.0;
        for (int n = 0; n < tc.poisNgbC; n++) {
            betas += datafn(tc.poisNgb[n].data()) * tc.poisCoef[n]; // * (tc.poisNgbDir ? 1.0 : -1.0); // Todo: Rho?

            if (cl.id().id == target) {
                std::cout << "\tNgb N=" << n << " " << tc.poisNgb[n]
                          << ", P=" << datafn(tc.poisNgb[n].data())
                          << ", b=" << tc.poisCoef[n] << "\n";
            }
        }

        if (cl.id().id == target) {
            std::cout << "\tsum(betas) = " << betas << "\n";
        }

        // alphas = 0.0;
        // for (int d = 0; d < 2*DIM; d++) {
        //     alphas += tc.poisAlpha[d]; // Todo: Rho?
        //     std::cout << "\talpha[" << d << "] = " << tc.poisAlpha[d] << "\n";
        // }
        alphas = cl->cenCoef;

        temp = (betas - cellfn(cl)) / alphas;

        // if (tc.poisNgbC == 2)
        // datafn(cl.data()) = temp;
        // else
        datafn(cl.data()) = omega*temp + (1.0-omega)*datafn(cl.data());

        if (temp != temp) {
            throw;
        }

        if (cl.id().id == target) {
            std::cout << "\tsum(alphas) = " << alphas << "\n";
            std::cout << "\tf = " << cellfn(cl) << "\n";
            std::cout << "\tP' = " << temp << "\n";
            std::cout << "\tP = " << datafn(cl.data()) << "\n";
        }
    }
}


double resid(fnDataRef datafn, fnCell cellfn)
{
    TreeCell* tc;

    double max = 0.0;
    cell_t maxat;

    for (auto& cl : leaves) {
        tc = &cl.group->cells[cl.index];

        cl->res = -cellfn(cl) - cl->cenCoef * datafn(cl.data());

        for (int n = 0; n < tc->poisNgbC; n++) {
            cl->res += datafn(tc->poisNgb[n].data()) * tc->poisCoef[n];
        }

        if (std::fabs(cl->res) > std::fabs(max)) {
            max = cl->res;
            maxat = cl;
        }
    }

    std::cerr << "Max Residual = " << max << " @ " << maxat << "\n";

    return max;
}


} // namespace tftt
