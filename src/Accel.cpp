/**
 *  @file   Accel.cpp
 *  @brief  Accel's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-30
 ***********************************************/

#include <iostream>
#include "Accel.h"


Matrix Accel(double x, const Matrix& Y) {
    std::cout << "[DEBUG] Accel inicio" << std::endl;
    std::cout << "Y: " << Y.getFilas() << "x" << Y.getColumnas() << std::endl;
    std::cout << "eopdata: " << eopdata.getFilas() << "x" << eopdata.getColumnas() << std::endl;
    std::cout << "PC: " << PC.getFilas() << "x" << PC.getColumnas() << std::endl;

    std::cout << "[DEBUG] IERS..." << std::endl;
    IERSResult ier = IERS(eopdata, AuxParam.Mjd_UTC + x/86400.0, 'l');
    std::cout << "[DEBUG] timediff..." << std::endl;
    TimeDiffResult td = timediff(ier.UT1_UTC, ier.TAI_UTC);

    double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400.0 + ier.UT1_UTC/86400.0;
    double Mjd_TT  = AuxParam.Mjd_UTC + x/86400.0 + td.TT_UTC/86400.0;

    std::cout << "[DEBUG] PrecMatrix..." << std::endl;
    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
    std::cout << "[DEBUG] NutMatrix..." << std::endl;
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    std::cout << "[DEBUG] PoleMatrix..." << std::endl;
    Matrix E = PoleMatrix(ier.x_pole, ier.y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    std::cout << "[DEBUG] JPL_Eph_DE430..." << std::endl;
    Ephemeris eph = JPL_Eph_DE430(MJD_TDB);

    std::cout << "[DEBUG] Construyendo r_sat..." << std::endl;
    Matrix r_sat(3,1);
    for (int i = 1; i <= 3; ++i)
        r_sat(i, 1) = Y(i, 1);

    std::cout << "[DEBUG] AccelHarmonic..." << std::endl;
    Matrix a = AccelHarmonic(r_sat, E, AuxParam.n, AuxParam.m);

    if (AuxParam.sun) {
        std::cout << "[DEBUG] AccelPointMass Sun..." << std::endl;
        a = a + AccelPointMass(r_sat, eph.r_Sun, GM_Sun);
    }
    if (AuxParam.moon) {
        std::cout << "[DEBUG] AccelPointMass Moon..." << std::endl;
        a = a + AccelPointMass(r_sat, eph.r_Moon, GM_Moon);
    }

    if (AuxParam.planets) {
        std::cout << "[DEBUG] AccelPointMass Planets..." << std::endl;
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

    std::cout << "[DEBUG] Accel fin OK" << std::endl;
    return dY;
}
