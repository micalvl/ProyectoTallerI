/**
 *  @file   NutMatrix.h
 *  @brief  NutMatrix method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/
#ifndef PROYECTOTALLERI_NUTMATRIX_H
#define PROYECTOTALLERI_NUTMATRIX_H

#include "Matrix.h"


/**
 * @brief Transformation from mean to true equator and equinox
 * @param[in] Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Nutation matrix
 */
Matrix NutMatrix(double Mjd_TT);


#endif //PROYECTOTALLERI_NUTMATRIX_H
