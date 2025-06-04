/**
 *  @file   gibbs.cpp
 *  @brief  gibbs method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-14
 ***********************************************/

#include "../include/gibbs.h"
#include "Sat_const.h"
#include "angl.h"
#include "unit.h"
#include <cmath>
using namespace std;




GibbsResult gibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3) {

    const double small = 1e-8;
    GibbsResult res;
    res.v2    = Matrix::zeros(3,1);
    res.theta = res.theta1 = res.copa = 0.0;
    res.error = "ok";

    double magr1 = r1.norm();
    double magr2 = r2.norm();
    double magr3 = r3.norm();

    Matrix p = Matrix::cross(r2, r3);
    Matrix q = Matrix::cross(r3, r1);
    Matrix w = Matrix::cross(r1, r2);

    Matrix pn  = unit(p);
    Matrix r1n = unit(r1);
    res.copa = asin(Matrix::dot(pn, r1n));
    if (fabs(Matrix::dot(r1n, pn)) > 0.017452406) {
        res.error = "not coplanar";
        return res;
    }

    Matrix d = p + q + w;
    double magd = d.norm();
    Matrix n = p.opsc(magr1) + q.opsc(magr2) + w.opsc(magr3);
    double magn = n.norm();
    Matrix nn = unit(n);
    Matrix dn = unit(d);

    if (magd < small || magn < small || Matrix::dot(nn, dn) < small) {
        res.error = "impossible";
        return res;
    }

    res.theta  = angl(r1, r2);
    res.theta1 = angl(r2, r3);

    double r1mr2 = magr1 - magr2;
    double r3mr1 = magr3 - magr1;
    double r2mr3 = magr2 - magr3;
    Matrix s = r3.opsc(r1mr2)
               + r2.opsc(r3mr1)
               + r1.opsc(r2mr3);
    Matrix b = Matrix::cross(d, r2);
    double l = sqrt(GM_Earth / (magd * magn));
    double tover2 = l / magr2;
    res.v2 = b.opsc(tover2) + s.opsc(l);

    return res;
}
