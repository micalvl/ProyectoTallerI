/**
 *  @file   IERS.h
 *  @brief  IERS method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-13
 ***********************************************/

#ifndef PROYECTOTALLERI_IERS_H
#define PROYECTOTALLERI_IERS_H


#include "Matrix.h"
#include "Sat_const.h"
#include <cmath>

struct IERSResult {
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;
};


/**
 * @brief Management of IERS time and polar motion data.
 * @param[in] eop EOP table
 * @param[in] Mjd_UTC Modified Julian Date (UTC) at which to retrieve/interpolate values.
 * @param[in] interp If 'l', perform linear interpolation between day floor(Mjd_UTC).
 * @return r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
 */
IERSResult IERS(const Matrix& eop, double Mjd_UTC, char interp = 'n');

#endif //PROYECTOTALLERI_IERS_H
