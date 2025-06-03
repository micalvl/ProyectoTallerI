//
// Created by micalvl on 03/04/2025.
//

#include "anglesdr.h"
#include "global.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Geodetic.h"
#include "LTC.h"
#include "doubler.h"
#include <cmath>
#include <string>
#include <iostream>

using namespace std;

AnglesDRResult anglesdr(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        const Matrix& rsite1,
        const Matrix& rsite2,
        const Matrix& rsite3
) {
    double magr1in = 1.1 * R_Earth;
    double magr2in = 1.11 * R_Earth;
    bool direct = true;
    double tol = 1e-8 * R_Earth;
    double pctchg = 0.005;

    double t1 = (Mjd1 - Mjd2) * 86400.0;
    double t3 = (Mjd3 - Mjd2) * 86400.0;

    Matrix los1(3,1), los2(3,1), los3(3,1);
    los1(1,1) = cos(el1) * sin(az1);
    los1(2,1) = cos(el1) * cos(az1);
    los1(3,1) = sin(el1);
    los2(1,1) = cos(el2) * sin(az2);
    los2(2,1) = cos(el2) * cos(az2);
    los2(3,1) = sin(el2);
    los3(1,1) = cos(el3) * sin(az3);
    los3(2,1) = cos(el3) * cos(az3);
    los3(3,1) = sin(el3);


    double lon1, lat1, h1;
    Geodetic(rsite1, lon1, lat1, h1);
    double lon2, lat2, h2;
    Geodetic(rsite2, lon2, lat2, h2);
    double lon3, lat3, h3;
    Geodetic(rsite3, lon3, lat3, h3);

    Matrix M1 = LTC(lon1, lat1);
    Matrix M2 = LTC(lon2, lat2);
    Matrix M3 = LTC(lon3, lat3);
    los1 = M1.transpose() * los1;
    los2 = M2.transpose() * los2;
    los3 = M3.transpose() * los3;

    auto transform = [](double Mjd_UTC, Matrix& los, Matrix& rsite) {
        IERSResult ier = IERS(eopdata, Mjd_UTC, 'l');
        TimeDiffResult td = timediff(ier.UT1_UTC, ier.TAI_UTC);
        double Mjd_TT = Mjd_UTC + td.TT_UTC / 86400.0;
        double Mjd_UT1 = Mjd_TT + (ier.UT1_UTC - td.TT_UTC) / 86400.0;
        Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
        Matrix N = NutMatrix(Mjd_TT);
        Matrix T = N * P;
        Matrix E = PoleMatrix(ier.x_pole, ier.y_pole) * GHAMatrix(Mjd_UT1) * T;
        los = E.transpose() * los;
        rsite = E.transpose() * rsite;
    };

    Matrix rsite1_t = rsite1;
    Matrix rsite2_t = rsite2;
    Matrix rsite3_t = rsite3;

    transform(Mjd1, los1, rsite1_t);
    transform(Mjd2, los2, rsite2_t);
    transform(Mjd3, los3, rsite3_t);

    double magr1old = 1e20, magr2old = 1e20;
    double magrsite1 = rsite1_t.norm(), magrsite2 = rsite2_t.norm(), magrsite3 = rsite3_t.norm();
    double cc1 = 2.0 * (los1.transpose() * rsite1_t)(1,1);
    double cc2 = 2.0 * (los2.transpose() * rsite2_t)(1,1);

    Matrix r2(3,1), r3(3,1), v2(3,1);
    double f1, f2, q1, magr1, magr2, a, deltae32;

    while (std::abs(magr1in - magr1old) > tol || std::abs(magr2in - magr2old) > tol) {
        auto res = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                           los1, los2, los3, rsite1_t, rsite2_t, rsite3_t, t1, t3, direct);
        r2 = res.r2; r3 = res.r3; f1 = res.f1; f2 = res.f2;
        magr1 = res.magr1; magr2 = res.magr2; a = res.a; deltae32 = res.deltae32;

        double f = 1.0 - a/magr2*(1.0 - cos(deltae32));
        double g = t3 - sqrt(a*a*a/GM_Earth)*(deltae32 - sin(deltae32));
        v2 = (r3 - r2.opsc(f)).divsc(g);

        double magr1o = magr1in;
        magr1in *= (1.0 + pctchg);
        double deltar1 = pctchg * magr1in;
        auto res1 = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2, los3, rsite1_t, rsite2_t, rsite3_t, t1, t3, direct);
        double pf1pr1 = (res1.f1 - f1) / deltar1, pf2pr1 = (res1.f2 - f2) / deltar1;

        magr1in = magr1o;
        double magr2o = magr2in;
        magr2in *= (1.0 + pctchg);
        double deltar2 = pctchg * magr2in;
        auto res2 = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2, los3, rsite1_t, rsite2_t, rsite3_t, t1, t3, direct);
        double pf1pr2 = (res2.f1 - f1) / deltar2, pf2pr2 = (res2.f2 - f2) / deltar2;

        double delta = pf1pr1 * pf2pr2 - pf2pr1 * pf1pr2;
        double delta1 = pf2pr2 * f1 - pf1pr2 * f2;
        double delta2 = pf1pr1 * f2 - pf2pr1 * f1;

        magr1old = magr1in; magr2old = magr2in;
        magr1in += -delta1 / delta;
        magr2in += -delta2 / delta;
    }

    auto res_final = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                             los1, los2, los3, rsite1_t, rsite2_t, rsite3_t, t1, t3, direct);
    r2 = res_final.r2; r3 = res_final.r3; a = res_final.a; deltae32 = res_final.deltae32;
    double f = 1.0 - a / magr2 * (1.0 - cos(deltae32));
    double g = t3 - sqrt(a*a*a / GM_Earth) * (deltae32 - sin(deltae32));
    v2 = (r3 - r2.opsc(f)).divsc(g);

    return {r2, v2};
}


