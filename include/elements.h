/**
 *  @file   elements.h
 *  @brief  elements method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_ELEMENTS_H
#define PROYECTOTALLERI_ELEMENTS_H


#include "Matrix.h"
#include "Sat_const.h"
#include "unit.h"
#include "angl.h"
#include <cmath>

struct ElementsResult {
    double p;
    double a;
    double e;
    double i;
    double Omega;
    double omega;
    double M;
};


/**
 * @brief Computes the osculating Keplerian elements from the satellite state vector for elliptic orbits.
 * @param[in] y State vector (x,y,z,vx,vy,vz)
 * @return       ElementsResult containing:
 *                - p     Semilatus rectum [m]
 *                - a     Semi-major axis
 *                - e     Eccentricity
 *                - i     Inclination [rad]
 *                - Omega Longitude of ascending node [rad]
 *                - omega Argument of pericenter [rad]
 *                - M     Mean anomaly [rad]
 */
ElementsResult elements(const Matrix& y);

#endif //PROYECTOTALLERI_ELEMENTS_H
