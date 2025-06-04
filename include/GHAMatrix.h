/**
 *  @file   GHAMatrix.h
 *  @brief  GHAMatrix method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_GHAMATRIX_H
#define PROYECTOTALLERI_GHAMATRIX_H

#include "Matrix.h"


/**
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system.
 * @param[in] Mjd_UT1 Modified Julian Date UT1.
 * @return Greenwich Hour Angle matrix.
 */
Matrix GHAMatrix(double Mjd_UT1);


#endif //PROYECTOTALLERI_GHAMATRIX_H
