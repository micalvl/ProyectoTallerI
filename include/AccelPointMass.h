/**
 *  @file   AccelPointMass.h
 *  @brief  AccelPointMass's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-13
 ***********************************************/

#include <array>
#include <cmath>
#include "../include/Matrix.h"

using namespace std;

#ifndef PROYECTOTALLERI_ACCELPOINTMASS_H
#define PROYECTOTALLERI_ACCELPOINTMASS_H


/**
 * @brief Computes the perturbational acceleration due to a point mass.
 * @param[in] r Satellite position vector.
 * @param[in] s Point mass position vector.
 * @param[in] GM Gravitational coefficient of point mass.
 * @return Acceleration (a=d^2r/dt^2).
 */
Matrix AccelPointMass(const Matrix& r, const Matrix& s, double GM);


#endif //PROYECTOTALLERI_ACCELPOINTMASS_H
