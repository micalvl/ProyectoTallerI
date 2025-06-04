/**
 *  @file   Accel.h
 *  @brief  Accel's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-30
 ***********************************************/

#ifndef PROYECTOTALLERI_ACCEL_H
#define PROYECTOTALLERI_ACCEL_H


#include "global.h"
#include "Matrix.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Mjday_TDB.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "Sat_const.h"

using namespace std;


/**
 * @brief Computes the acceleration of an Earth orbiting satellite due to
 *   - the Earth's harmonic gravity field,
 *   - the gravitational perturbations of the Sun and Moon
 *   - the solar radiation pressure and
 *   - the atmospheric drag
 * @param[in] x Terrestrial Time (Modified Julian Date).
 * @param[in] Y Satellite state vector in the ICRF/EME2000 system.
 * @return Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system.
 */
Matrix Accel(double x, const Matrix& Y);

#endif //PROYECTOTALLERI_ACCEL_H
