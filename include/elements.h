//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_ELEMENTS_H
#define PROYECTOTALLERI_ELEMENTS_H


#include "Matrix.h"
#include "Sat_const.h"
#include "unit.h"
#include "angl.h"
#include <cmath>

struct ElementsResult {
    double p;      ///< semilatus rectum [m]
    double a;      ///< semi-major axis [m]
    double e;      ///< eccentricity
    double i;      ///< inclination [rad]
    double Omega;  ///< longitude of ascending node [rad]
    double omega;  ///< argument of pericenter [rad]
    double M;      ///< mean anomaly [rad]
};

ElementsResult elements(const Matrix& y);

#endif //PROYECTOTALLERI_ELEMENTS_H
