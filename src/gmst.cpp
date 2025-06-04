/**
 *  @file   gmst.cpp
 *  @brief  gmst function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-23
 ***********************************************/

#include <cmath>
#include "../include/gmst.h"
#include "../include/Sat_const.h"
#include "Frac.h"


double gmst(double Mjd_UT1){
    double gmstime, Secs, MJD_J2000, Mjd_0, UT1, T, T_0, gmst;
    Secs = 86400.0;                       // Seconds per day
    MJD_J2000= 51544.5;

    Mjd_0 = floor(Mjd_UT1);
    UT1   = Secs*(Mjd_UT1-Mjd_0);         // [s]
    T_0   = (Mjd_0  -MJD_J2000)/36525.0;
    T     = (Mjd_UT1-MJD_J2000)/36525.0;

    gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1
    + (0.093104-6.2e-6*T)*T*T;    // [s]

    gmstime = 2*M_PI*Frac(gmst/Secs);      // [rad], 0..2pi
    return gmstime;
}
