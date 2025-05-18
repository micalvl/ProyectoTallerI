//
// Created by micalvl on 03/04/2025.
//

#include "../include/hgibbs.h"

#include <cmath>
using namespace std;

HGibbsResult hgibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3,
                    double Mjd1, double Mjd2, double Mjd3)
{
    HGibbsResult res;
    res.v2 = Matrix(3,1);
    res.theta = res.theta1 = res.copa = 0.0;
    res.error = "ok";

    double magr1 = r1.norm();
    double magr2 = r2.norm();
    double magr3 = r3.norm();
    const double tolangle = 0.01745329251994;

    double dt21 = (Mjd2 - Mjd1) * 86400.0;
    double dt31 = (Mjd3 - Mjd1) * 86400.0;
    double dt32 = (Mjd3 - Mjd2) * 86400.0;

    Matrix p = Matrix::cross(r2, r3);
    Matrix pn = unit(p);
    Matrix r1n = unit(r1);
    res.copa = asin(Matrix::dot(pn,r1n));

    if (fabs(Matrix::dot(r1n, pn)) > 0.017452406) {
        res.error = "not coplanar";
        return res;
    }

    res.theta  = angl(r1, r2);
    res.theta1 = angl(r2, r3);


    if (res.theta > tolangle || res.theta1 > tolangle) {
        res.error = "angl > 1ø";
        return res;
    }

    double term1 = -dt32 * (1.0/(dt21*dt31) +
                            GM_Earth/(12.0 * magr1*magr1*magr1));
    double term2 = (dt32 - dt21) * (1.0/(dt21*dt32) +
                                    GM_Earth/(12.0 * magr2*magr2*magr2));
    double term3 =  dt21 * (1.0/(dt32*dt31) +
                            GM_Earth/(12.0 * magr3*magr3*magr3));

    res.v2 = r1.opsc(term1)
             + r2.opsc(term2)
             + r3.opsc(term3);

    return res;
}