//
// Created by micalvl on 03/04/2025.
//
#include "anglesg.h"
#include "global.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Geodetic.h"
#include "LTC.h"
#include "gibbs.h"
#include "hgibbs.h"
#include "elements.h"
#include "angl.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>
#include <complex>
#include <iostream>
#include "../Rpoly/Rpoly/rpoly.h"

AnglesGResult anglesg(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        const Matrix& Rs1, const Matrix& Rs2, const Matrix& Rs3
) {
    Matrix L1(3,1), L2(3,1), L3(3,1);
    L1(1,1) = cos(el1)*sin(az1); L1(2,1) = cos(el1)*cos(az1); L1(3,1) = sin(el1);
    L2(1,1) = cos(el2)*sin(az2); L2(2,1) = cos(el2)*cos(az2); L2(3,1) = sin(el2);
    L3(1,1) = cos(el3)*sin(az3); L3(2,1) = cos(el3)*cos(az3); L3(3,1) = sin(el3);

    double lon1, lat1, h1, lon2, lat2, h2, lon3, lat3, h3;
    Geodetic(Rs1, lon1, lat1, h1);
    Geodetic(Rs2, lon2, lat2, h2);
    Geodetic(Rs3, lon3, lat3, h3);

    Matrix M1 = LTC(lon1, lat1);
    Matrix M2 = LTC(lon2, lat2);
    Matrix M3 = LTC(lon3, lat3);

    Matrix Lb1 = M1.transpose().operator*(L1);

    Matrix Lb2 = M2.transpose().operator*(L2);

    Matrix Lb3 = M3.transpose().operator*(L3);

    auto transform = [](double Mjd_UTC, Matrix& L, Matrix& Rs) {
        IERSResult ier = IERS(eopdata, Mjd_UTC, 'l');
        TimeDiffResult td = timediff(ier.UT1_UTC, ier.TAI_UTC);
        double Mjd_TT = Mjd_UTC + td.TT_UTC/86400.0;
        double Mjd_UT1 = Mjd_TT + (ier.UT1_UTC - td.TT_UTC)/86400.0;
        Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
        Matrix N = NutMatrix(Mjd_TT);
        Matrix T = N * P;
        Matrix E = PoleMatrix(ier.x_pole, ier.y_pole) * GHAMatrix(Mjd_UT1) * T;
        L = E.transpose() * L;
        Rs = E.transpose() * Rs;
    };

    Matrix Rs1_t = Rs1;
    Matrix Rs2_t = Rs2;
    Matrix Rs3_t = Rs3;

    transform(Mjd1, Lb1, Rs1_t);
    transform(Mjd2, Lb2, Rs2_t);
    transform(Mjd3, Lb3, Rs3_t);

    double tau1 = (Mjd1 - Mjd2) * 86400.0;
    double tau3 = (Mjd3 - Mjd2) * 86400.0;
    double a1 = tau3 / (tau3 - tau1);
    double a3 = -tau1 / (tau3 - tau1);
    double b1 = tau3 / (6*(tau3 - tau1)) * ((tau3 - tau1)*(tau3 - tau1) - tau3*tau3);
    double b3 = -tau1 / (6*(tau3 - tau1)) * ((tau3 - tau1)*(tau3 - tau1) - tau1*tau1);

    Matrix L(3,3), R(3,3);
    for (int i = 1; i <= 3; ++i) {
        L(i,1) = Lb1(i,1);
        L(i,2) = Lb2(i,1);
        L(i,3) = Lb3(i,1);
    }
    for (int i = 1; i <= 3; ++i) {
        R(i,1) = Rs1_t(i,1);
        R(i,2) = Rs2_t(i,1);
        R(i,3) = Rs3_t(i,1);
    }
    L.print();
    R.print();

    Matrix Linv = L.inverse();
    Linv.print();

    Matrix D = Linv * R;
    D.print();
    double Ccye = 2.0 * (Lb2.transpose() * Rs2_t)(1,1);
    double d1s = D(2,1)*a1 - D(2,2) + D(2,3)*a3;
    double d2s = D(2,1)*b1 + D(2,3)*b3;

    vector<double> poly(9, 0.0);
    poly[0] = 1.0;
    poly[2] = -(d1s*d1s + d1s*Ccye + Rs2_t.norm()*Rs2_t.norm());
    poly[6] = -GM_Earth*(d2s*Ccye + 2*d1s*d2s);
    poly[8] = -GM_Earth*GM_Earth*d2s*d2s;

    // RPOLY
    double coef[9];
    for (int i = 0; i < 9; ++i) coef[i] = poly[i];

    double zeror[8], zeroi[8];
    int nroots = real_poly_roots(coef, 8, zeror, zeroi);

    double bigr2 = -numeric_limits<double>::max();
    for (int i = 0; i < nroots; ++i) {
        if (fabs(zeroi[i]) < 1e-8 && zeror[i] > bigr2) {
            bigr2 = zeror[i];
        }
    }

    if (bigr2 <= 0) throw runtime_error("anglesg: Error");

    double u = GM_Earth / (bigr2 * bigr2 * bigr2);
    double C1 = a1 + b1 * u, C2 = -1, C3 = a3 + b3 * u;
    Matrix C(3, 1);
    C(1, 1) = C1; C(2, 1) = C2; C(3, 1) = C3;

    Matrix temp = D * C;
    temp = temp.opsc(-1.0);

    double rho1 = temp(1,1)/(a1 + b1*u);
    double rho2 = -temp(2,1);
    double rho3 = temp(3,1)/(a3 + b3*u);

    Matrix r1 = Rs1_t + Lb1.opsc(rho1);
    Matrix r2 = Rs2_t + Lb2.opsc(rho2);
    Matrix r3 = Rs3_t + Lb3.opsc(rho3);

    Matrix v2(3,1);
    std::string error;
    GibbsResult res = gibbs(r1, r2, r3);
    v2 = res.v2;
    error = res.error;
    if (error != "ok") {
        HGibbsResult res2 = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);
        v2 = res2.v2;
        error = res2.error;
    }

    return {r2, v2};
}