/**
 *  @file   NutAngles.h
 *  @brief  NutAngles function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-23
 ***********************************************/

#ifndef PROYECTOTALLERI_NUTANGLES_H
#define PROYECTOTALLERI_NUTANGLES_H


#include "Matrix.h"
#include "Sat_const.h"
#include <cmath>


/**
 * @brief Nutation in longitude and obliquity
 * @param[in] Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return dpsi,deps  Nutation Angles
 */
Matrix NutAngles(double Mjd_TT);


#endif //PROYECTOTALLERI_NUTANGLES_H
