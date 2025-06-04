/**
 *  @file   Geodetic.h
 *  @brief  Geodetic method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/


#ifndef PROYECTOTALLERI_GEODETIC_H
#define PROYECTOTALLERI_GEODETIC_H


#include "Matrix.h"


/**
 * @brief geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from given position vector (r [m])
 * @param[in] r Matrix
 * @param[out] lon Longitude [rad]
 * @param[out] lat Latitude [rad]
 * @param[out] h Altitude [m]
 */
void Geodetic(const Matrix& r, double& lon, double& lat, double& h);


#endif //PROYECTOTALLERI_GEODETIC_H
