/**
 *  @file   angl.h
 *  @brief  angl method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-13
 ***********************************************/

#ifndef PROYECTOTALLERI_ANGL_H
#define PROYECTOTALLERI_ANGL_H


#include "Matrix.h"
#include <cmath>
using namespace std;


/**
 * @brief Computes the angle between two vectors.
 * @param[in] vec1  First input vector.
 * @param[in] vec2  Second input vector.
 * @return (double) angle between the two vectors  -pi to pi.
 */
double angl(const Matrix& vec1, const Matrix& vec2);


#endif //PROYECTOTALLERI_ANGL_H
