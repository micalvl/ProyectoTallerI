//
// Created by micalvl on 03/04/2025.
//

#include "doubler.h"
#include <cmath>
using namespace std;

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


    double D1 = cc1*cc1 - 4*(magrsite1*magrsite1 - magr1in*magr1in);
    double D2 = cc2*cc2 - 4*(magrsite2*magrsite2 - magr2in*magr2in);
    double rho1 = (-cc1 + std::sqrt(D1))/2.0;
    double rho2 = (-cc2 + std::sqrt(D2))/2.0;

    Matrix r1 = los1.opsc(rho1) + rsite1;
    R.r2       = los2.opsc(rho2) + rsite2;
    R.magr1    = r1.norm();
    R.magr2    = R.r2.norm();


    Matrix w = Matrix::cross(r1, R.r2).opsc(1.0/(R.magr1*R.magr2));
    if (!direct) w = w.opsc(-1.0);


    double rho3 = -Matrix::dot(rsite3, w) / Matrix::dot(los3, w);
    R.r3 = los3.opsc(rho3) + rsite3;
    double magr3 = R.r3.norm();


    auto angle = [&](const Matrix& A, const Matrix& B, double mA, double mB){
        double c = Matrix::dot(A,B)/(mA*mB);
        double s = Matrix::cross(A,B).norm()/(mA*mB);
        return atan2(s, c);
    };
    double dv21 = angle(R.r2, r1,    R.magr2, R.magr1);
    double dv31 = angle(R.r3, r1,    magr3,   R.magr1);
    double dv32 = angle(R.r3, R.r2,  magr3,   R.magr2);


    double c1,c3,p;
    if (dv31 > M_PI) {
        c1 = (R.magr2*sin(dv32))/(R.magr1*sin(dv31));
        c3 = (R.magr2*sin(dv21))/(magr3*sin(dv31));
        p  = (c1*R.magr1 + c3*magr3 - R.magr2)/(c1 + c3 - 1);
    } else {
        c1 = (R.magr1*sin(dv31))/(R.magr2*sin(dv32));
        c3 = (R.magr1*sin(dv21))/(magr3*std::sin(dv32));
        p  = (c3*magr3 - c1*R.magr2 + R.magr1)/(-c1 + c3 + 1);
    }


    double ecosv2 = p/R.magr2 - 1;
    double esinv2 = (-cos(dv21)*ecosv2 + (p/R.magr1 - 1))/sin(dv21);
    double e = sqrt(ecosv2*ecosv2 + esinv2*esinv2);
    R.a = p/(1 - e*e);


    double n = sqrt(GM_Earth/(R.a*R.a*R.a));
    double s = R.magr2/p * sqrt(1-e*e) * esinv2;
    double c = R.magr2/p * (e*e + ecosv2);

    double sinde32 = magr3/sqrt(R.a*p)*std::sin(dv32) - magr3/p*(1-cos(dv32))*s;
    R.deltae32 = atan2(sinde32, 1 - R.magr2*magr3/(R.a*p)*(1-cos(dv32)));

    double sinde21 = R.magr1/sqrt(R.a*p)*std::sin(dv21) + R.magr1/p*(1-cos(dv21))*s;
    double deltae21= atan2(sinde21, 1 - R.magr2*R.magr1/(R.a*p)*(1-cos(dv21)));

    double deltam32 = R.deltae32 + 2*s*pow(sin(R.deltae32/2),2) - c*sin(R.deltae32);
    double deltam12 = -deltae21 + 2*s*pow(sin(deltae21/2),2) + c*sin(deltae21);

    R.f1 = t1 - deltam12/n;
    R.f2 = t3 - deltam32/n;
    R.q1 = sqrt(R.f1*R.f1 + R.f2*R.f2);

    return R;
}
