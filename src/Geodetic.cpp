//
// Created by micalvl on 03/04/2025.
//

#include "../include/Geodetic.h"
#include "Sat_const.h"
#include <cmath>
using namespace std;

void Geodetic(const Matrix& r, double& lon, double& lat, double& h) {
    const double R_equ   = R_Earth;
    const double f       = f_Earth;
    const double epsRequ = 1e-12 * R_equ;
    const double e2      = f * (2.0 - f);

    double X    = r(1,1);
    double Y    = r(2,1);
    double Z    = r(3,1);
    double rho2 = X*X + Y*Y;

    double dZ   = e2 * Z;
    double ZdZ, Nh, SinPhi, N, dZ_new;
    while (true) {
        ZdZ    = Z + dZ;
        Nh     = sqrt(rho2 + ZdZ*ZdZ);
        SinPhi = ZdZ / Nh;
        N      = R_equ / sqrt(1.0 - e2 * SinPhi * SinPhi);
        dZ_new = N * e2 * SinPhi;
        if (std::fabs(dZ - dZ_new) < epsRequ)
            break;
        dZ = dZ_new;
    }

    lon = atan2(Y, X);
    lat = atan2(ZdZ, sqrt(rho2));
    h   = Nh - N;
}