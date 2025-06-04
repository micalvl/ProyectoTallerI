/**
 *  @file   Position.h
 *  @brief  Position method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-20
 ***********************************************/

#ifndef PROYECTOTALLERI_POSITION_H
#define PROYECTOTALLERI_POSITION_H

#include "Matrix.h"
#include "Sat_const.h"
#include <cmath>
using namespace std;


/**
 * @brief Computes the geocentric position vector from geodetic coordinates.
 * @param[in] lon Longitude [rad].
 * @param[in] lat Latitude [rad].
 * @param[in] h Altitude [m].
 * @return 3×1 position vector.
 */
Matrix Position(double lon, double lat, double h);


#endif //PROYECTOTALLERI_POSITION_H
