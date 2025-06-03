/**
 *  @file   G_AccelHarmonic.h
 *  @brief  G_AccelHarmonic function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-07
 ***********************************************/

#ifndef PROYECTOTALLERI_G_ACCELHARMONIC_H
#define PROYECTOTALLERI_G_ACCELHARMONIC_H


#include "Matrix.h"
#include "AccelHarmonic.h"
/**
 * @brief Computes the gradient of Earth's harmonic gravitational field.
 *
 * @param r Position vector of the satellite in the true-of-date reference frame.
 * @param U Transformation matrix to the Earth-fixed frame.
 * @param n_max Maximum degree of the gravitational model.
 * @param m_max Maximum order of the gravitational model.
 * @return 3x3 matrix with the gradient (G = da/dr) in the true-of-date reference frame.
 */
Matrix G_AccelHarmonic(const Matrix& r, const Matrix& E, int n_max, int m_max);


#endif //PROYECTOTALLERI_G_ACCELHARMONIC_H
