/**
 *  @file   elements.cpp
 *  @brief  elements method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#include "../include/elements.h"

ElementsResult elements(const Matrix& y) {

    Matrix r(3,1), v(3,1);
    for (int k = 1; k <= 3; ++k) {
        r(k,1) = y(k,1);
        v(k,1) = y(k+3,1);
    }

    Matrix h = Matrix::cross(r, v);
    double magh = h.norm();
    double p = magh * magh / GM_Earth;


    double Omega = atan2(h(1,1), -h(2,1));
    Omega = fmod(Omega, 2*M_PI);
    if (Omega < 0) Omega += 2*M_PI;

    double i = atan2(sqrt(h(1,1)*h(1,1) + h(2,1)*h(2,1)), h(3,1));

    double u = atan2(r(3,1) * magh, -r(1,1)*h(2,1) + r(2,1)*h(1,1));

    double R  = r.norm();
    double dotv2 = Matrix::dot(v, v);

    double a = 1.0 / (2.0/R - dotv2/GM_Earth);

    double eCosE = 1.0 - R/a;
    double eSinE = Matrix::dot(r, v) / sqrt(GM_Earth * a);

    double e2 = eCosE*eCosE + eSinE*eSinE;
    double e  = sqrt(e2);
    double E  = atan2(eSinE, eCosE);

    double M = fmod(E - eSinE, 2*M_PI);
    if (M < 0) M += 2*M_PI;

    double nu = atan2(sqrt(1.0 - e2)*eSinE, eCosE - e2);

    double omega = fmod(u -nu, 2*M_PI);
    if (omega < 0) omega += 2*M_PI;

    return { p, a, e, i, Omega, omega, M };
}