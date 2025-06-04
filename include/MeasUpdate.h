/**
 *  @file   MeasUpdate.h
 *  @brief  MeasUpdate method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-20
 ***********************************************/

#ifndef PROYECTOTALLERI_MEASUPDATE_H
#define PROYECTOTALLERI_MEASUPDATE_H

#include "Matrix.h"


/**
 * @brief Updates state and covariance with a measurement.
 * @param[in,out] s Matrix
 * @param[in] z Matrix
 * @param[in] g Matrix
 * @param[in] s Matrix
 * @param[in] G Matrix
 * @param[in,out] P Matrix
 * @param[in] n dimension
 * @return MeasUpdate Matrix
 */
Matrix MeasUpdate(
        Matrix&       x,
        const Matrix& z,
        const Matrix& g,
        const Matrix& s,
        const Matrix& G,
        Matrix&       P,
        int           n
);


#endif //PROYECTOTALLERI_MEASUPDATE_H
