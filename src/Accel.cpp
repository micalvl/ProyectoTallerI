/**
 *  @file   Accel.cpp
 *  @brief  Accel's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-30
 ***********************************************/

#include "Accel.h"


Matrix Accel(double x, const Matrix& Y) {

    IERSResult ier = IERS(eopdata, AuxParam.Mjd_UTC + x/86400.0, 'l');
    double x_pole  = ier.x_pole;
    double y_pole  = ier.y_pole;
    double UT1_UTC = ier.UT1_UTC;
    double LOD     = ier.LOD;
    double dpsi    = ier.dpsi;
    double deps    = ier.deps;
    double dx_pole = ier.dx_pole;
    double dy_pole = ier.dy_pole;
    double TAI_UTC = ier.TAI_UTC;

    TimeDiffResult td = timediff(UT1_UTC, TAI_UTC);
    double UT1_TAI = td.UT1_TAI;
    double UTC_GPS = td.UTC_GPS;
    double UT1_GPS = td.UT1_GPS;
    double TT_UTC  = td.TT_UTC;
    double GPS_UTC = td.GPS_UTC;

    double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400.0 + UT1_UTC/86400.0;
    double Mjd_TT  = AuxParam.Mjd_UTC + x/86400.0 + TT_UTC/86400.0;

    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    Ephemeris eph = JPL_Eph_DE430(MJD_TDB);
    Matrix r_Mercury = eph.r_Mercury;
    Matrix r_Venus   = eph.r_Venus;
    Matrix r_Earth   = eph.r_Earth;
    Matrix r_Mars    = eph.r_Mars;
    Matrix r_Jupiter = eph.r_Jupiter;
    Matrix r_Saturn  = eph.r_Saturn;
    Matrix r_Uranus  = eph.r_Uranus;
    Matrix r_Neptune = eph.r_Neptune;
    Matrix r_Pluto   = eph.r_Pluto;
    Matrix r_Moon    = eph.r_Moon;
    Matrix r_Sun     = eph.r_Sun;

    Matrix r_sat(3,1);
    r_sat(1,1) = Y(1,1);
    r_sat(2,1) = Y(2,1);
    r_sat(3,1) = Y(3,1);

    Matrix a = AccelHarmonic(r_sat, E, AuxParam.n, AuxParam.m);

    if (AuxParam.sun) {
        a = a + AccelPointMass(r_sat, r_Sun,  GM_Sun);
    }
    if (AuxParam.moon) {
        a = a + AccelPointMass(r_sat, r_Moon, GM_Moon);
    }

    if (AuxParam.planets) {
        a = a + AccelPointMass(r_sat, r_Mercury, GM_Mercury);
        a = a + AccelPointMass(r_sat, r_Venus,   GM_Venus);
        a = a + AccelPointMass(r_sat, r_Mars,    GM_Mars);
        a = a + AccelPointMass(r_sat, r_Jupiter, GM_Jupiter);
        a = a + AccelPointMass(r_sat, r_Saturn,  GM_Saturn);
        a = a + AccelPointMass(r_sat, r_Uranus,  GM_Uranus);
        a = a + AccelPointMass(r_sat, r_Neptune, GM_Neptune);
        a = a + AccelPointMass(r_sat, r_Pluto,   GM_Pluto);
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
