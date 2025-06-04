/**
 *  @file   PoleMatrix.h
 *  @brief  PoleMatrix method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_POLEMATRIX_H
#define PROYECTOTALLERI_POLEMATRIX_H


#include "Matrix.h"


/**
 * @brief Constructs the pole matrix transforming from pseudo Earth-fixed to Earth-fixed.
 * @param[in] xp Polar motion coordinate x (rad).
 * @param[in] yp Polar motion coordinate y (rad).
 * @return Pole matrix
 */
Matrix PoleMatrix(double xp, double yp);


#endif //PROYECTOTALLERI_POLEMATRIX_H
