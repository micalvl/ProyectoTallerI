/**
 *  @file   PrecMatrix.h
 *  @brief  PrecMatrix method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_PRECMATRIX_H
#define PROYECTOTALLERI_PRECMATRIX_H


#include "Matrix.h"

/**
 * @brief Precession transformation of equatorial coordinates.
 * @param[in]  Mjd_1  Epoch given (Modified Julian Date TT).
 * @param[in]  Mjd_2  Epoch to precess to (Modified Julian Date TT).
 * @return Precession transformation matrix.
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif //PROYECTOTALLERI_PRECMATRIX_H
