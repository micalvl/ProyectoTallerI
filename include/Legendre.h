/**
 *  @file   Legendre.h
 *  @brief  Legendre method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-10
 ***********************************************/

#ifndef PROYECTOTALLERI_LEGENDRE_H
#define PROYECTOTALLERI_LEGENDRE_H


#include "Matrix.h"


/**
 * @brief Computes associated Legendre functions.
 * @param[in] n Maximum degree.
 * @param[in] m Maximum order.
 * @param[in] phi Latitude angle [rad].
 * @param[out] pnm pnm matrix.
 * @param[out] dpnm dpnm matrix.
 */
void Legendre(int n,int m, double fi, Matrix &pnm, Matrix &dpnm);


#endif //PROYECTOTALLERI_LEGENDRE_H
