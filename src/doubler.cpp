//
// Created by micalvl on 03/04/2025.
//

#include "../include/doubler.h"
#include "Sat_const.h"
#include <cmath>

DoublerResult doubler(
        double cc1, double cc2,
        double magrsite1, double magrsite2,
        double magr1in,   double magr2in,
        const Matrix& los1, const Matrix& los2, const Matrix& los3,
        const Matrix& rsite1,const Matrix& rsite2,const Matrix& rsite3,
        double t1, double t3,
        bool direct
) {
    DoublerResult R;

    double D1 = cc1*cc1 - 4.0*(magrsite1*magrsite1 - magr1in*magr1in);
    double D2 = cc2*cc2 - 4.0*(magrsite2*magrsite2 - magr2in*magr2in);
    double rho1 = (-cc1 + sqrt(D1))/2.0;
    double rho2 = (-cc2 + sqrt(D2))/2.0;

    Matrix r1 = los1.opsc(rho1) + rsite1;
    R.r2 = los2.opsc(rho2) + rsite2;

    R.magr1 = r1.norm();
    R.magr2 = R.r2.norm();

    Matrix w = Matrix::cross(r1, R.r2).opsc(1.0/(R.magr1*R.magr2));
    if (!direct) {
        w = w.opsc(-1.0);
    }


    double rho3 = - Matrix::dot(rsite3, w) / Matrix::dot(los3, w);
    R.r3 = los3.opsc(rho3) + rsite3;
    double magr3 = R.r3.norm();

    double cosdv21 = Matrix::dot(R.r2, r1)/(R.magr2*R.magr1);
    double sindv21 = Matrix::cross(R.r2, r1).norm()/(R.magr2*R.magr1);
    double dv21 = atan2(sindv21, cosdv21);

    double cosdv31 = Matrix::dot(R.r3, r1)/(magr3*R.magr1);
    double sindv31 = sqrt(1.0 - cosdv31*cosdv31);
    double dv31 = atan2(sindv31, cosdv31);

    double cosdv32 = Matrix::dot(R.r3, R.r2)/(magr3*R.magr2);
    double sindv32 = Matrix::cross(R.r3, R.r2).norm()/(magr3*R.magr2);

    double c1, c3, p;
    if (dv31 > M_PI) {
        c1 = (R.magr2*sindv32)/(R.magr1*sindv31);
        c3 = (R.magr2*sindv21)/(magr3*sindv31);
        p = (c1*R.magr1 + c3*magr3 - R.magr2)/(c1 + c3 - 1.0);
    } else {
        c1 = (R.magr1*sindv31)/(R.magr2*sindv32);
        c3 = (R.magr1*sindv21)/(magr3   *sindv32);
        p = (c3*magr3 - c1*R.magr2 + R.magr1)/(-c1 + c3 + 1.0);
    }

    double ecosv1 = p/R.magr1 - 1.0;
    double ecosv2 = p/R.magr2 - 1.0;
    double ecosv3 = p/magr3   - 1.0;

    double esinv2;
    if (fabs(dv21 - M_PI) > 1e-12) {
        esinv2 = (-cosdv21*ecosv2 + ecosv1)/sindv21;
    } else {
        esinv2 = (cosdv32*ecosv2 - ecosv3)/sindv31;
    }

    double e = sqrt(ecosv2*ecosv2 + esinv2*esinv2);
    R.a = p/(1.0 - e*e);

    double n, s, c, deltae32, deltae21, deltam32, deltam12;
    if (e*e < 0.99) {
        n = sqrt(GM_Earth/(R.a*R.a*R.a));
        s = R.magr2/p * sqrt(1.0 - e*e) * esinv2;
        c = R.magr2/p * (e*e + ecosv2);

        double sinde32 = magr3/sqrt(R.a*p)*sindv32 - magr3/p*(1.0 - cosdv32)*s;
        double cosde32 = 1.0 - R.magr2*magr3/(R.a*p)*(1.0 - cosdv32);
        deltae32 = atan2(sinde32, cosde32);

        double sinde21 = R.magr1/sqrt(R.a*p)*sindv21 + R.magr1/p*(1.0 - cosdv21)*s;
        double cosde21 = 1.0 - R.magr2*R.magr1/(R.a*p)*(1.0 - cosdv21);
        deltae21 = atan2(sinde21, cosde21);

        deltam32 = deltae32 + 2.0*s*pow(sin(deltae32/2.0),2) - c*sin(deltae32);
        deltam12 = -deltae21 + 2.0*s*pow(sin(deltae21/2.0),2) + c*sin(deltae21);
    } else {
        n = sqrt(GM_Earth/(-R.a*R.a*R.a));
        s = R.magr2/p * sqrt(e*e - 1.0) * esinv2;
        c = R.magr2/p * (e*e + ecosv2);

        double sindh32 = magr3/sqrt(-R.a*p)*sindv32 - magr3/p*(1.0 - cosdv32)*s;
        double deltah32 = log(sindh32 + sqrt(sindh32*sindh32+1.0));

        double sindh21 = R.magr1/sqrt(-R.a*p)*sindv21 + R.magr1/p*(1.0 - cosdv21)*s;
        double deltah21 = log(sindh21 + sqrt(sindh21*sindh21+1.0));

        deltam32 = -deltah32 + 2.0*s*pow(sinh(deltah32/2.0),2) + c*sinh(deltah32);
        deltam12 =  deltah21 + 2.0*s*pow(sinh(deltah21/2.0),2) - c*sinh(deltah21);
    }

    R.f1 = t1 - deltam12/n;
    R.f2 = t3 - deltam32/n;
    R.q1 = sqrt(R.f1*R.f1 + R.f2*R.f2);

    return R;
}
