/**
 *  @file   Accel.cpp
 *  @brief  Accel's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-26
 ***********************************************/

#include <iostream>
#include "Accel.h"


Matrix Accel(double x, const Matrix& Y) {

    IERSResult ier = IERS(eopdata, AuxParam.Mjd_UTC + x/86400.0, 'l');
    TimeDiffResult td = timediff(ier.UT1_UTC, ier.TAI_UTC);

    double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400.0 + ier.UT1_UTC/86400.0;
    double Mjd_TT  = AuxParam.Mjd_UTC + x/86400.0 + td.TT_UTC/86400.0;

    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(ier.x_pole, ier.y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    Ephemeris eph = JPL_Eph_DE430(MJD_TDB);

    Matrix r_sat(3,1);
    for (int i = 1; i <= 3; ++i)
        r_sat(i, 1) = Y(i, 1);

    Matrix a = AccelHarmonic(r_sat, E, AuxParam.n, AuxParam.m);

    if (AuxParam.sun) {
        a = a + AccelPointMass(r_sat, eph.r_Sun, GM_Sun);
    }
    if (AuxParam.moon) {
        a = a + AccelPointMass(r_sat, eph.r_Moon, GM_Moon);
    }

    if (AuxParam.planets) {
        a = a + AccelPointMass(r_sat, eph.r_Mercury, GM_Mercury);
        a = a + AccelPointMass(r_sat, eph.r_Venus,   GM_Venus);
        a = a + AccelPointMass(r_sat, eph.r_Mars,    GM_Mars);
        a = a + AccelPointMass(r_sat, eph.r_Jupiter, GM_Jupiter);
        a = a + AccelPointMass(r_sat, eph.r_Saturn,  GM_Saturn);
        a = a + AccelPointMass(r_sat, eph.r_Uranus,  GM_Uranus);
        a = a + AccelPointMass(r_sat, eph.r_Neptune, GM_Neptune);
        a = a + AccelPointMass(r_sat, eph.r_Pluto,   GM_Pluto);
    }


    Matrix dY(6,1);
    dY(1,1) = Y(4,1);
    dY(2,1) = Y(5,1);
    dY(3,1) = Y(6,1);
    dY(4,1) = a(1,1);
    dY(5,1) = a(2,1);
    dY(6,1) = a(3,1);

    return dY;
}
