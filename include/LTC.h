/**
 *  @file   LTC.h
 *  @brief  LTC method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_LTC_H
#define PROYECTOTALLERI_LTC_H

#include "Matrix.h"


/**
 * @brief Transformation from Greenwich meridian system to local tangent coordinates.
 * @param[in]  lon  Geodetic East longitude [rad].
 * @param[in]  lat  Geodetic latitude [rad].
 * @return Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system.
 */
Matrix LTC(double lon, double lat);


#endif //PROYECTOTALLERI_LTC_H
