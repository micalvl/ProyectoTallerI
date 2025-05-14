/**
 *  @file   AccelHarmonic.h
 *  @brief  AccelHarmonic function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-23
 ***********************************************/

#ifndef PROYECTOTALLERI_ACCELHARMONIC_H
#define PROYECTOTALLERI_ACCELHARMONIC_H


#include "Matrix.h"


/**
 * @brief Computes the acceleration due to the harmonic gravity field of a central body
 *
 * @param r Satellite position vector in the inertial reference frame (3x1)
 * @param E Rotation matrix from inertial frame to body-fixed frame (3x3)
 * @param n_max Maximum degree of the gravity field model
 * @param m_max Maximum order of the gravity field model (must satisfy m_max <= n_max)
 * @return Acceleration vector in the inertial frame (3x1)
 */
Matrix AccelHarmonic(const Matrix& r, const Matrix& E, int n_max, int m_max);


#endif //PROYECTOTALLERI_ACCELHARMONIC_H
