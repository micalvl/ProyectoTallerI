/**
 *  @file   unit.h
 *  @brief  unit method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-14
 ***********************************************/

#ifndef PROYECTOTALLERI_UNIT_H
#define PROYECTOTALLERI_UNIT_H

#include "Matrix.h"

/**
 * @brief Computes the unit (normalized) vector.
 * @param[in] vec Input vector (any size).
 * @return A vector of the same dimensions, with each element divided by the vector’s norm.
 */
Matrix unit(const Matrix& vec);


#endif //PROYECTOTALLERI_UNIT_H
