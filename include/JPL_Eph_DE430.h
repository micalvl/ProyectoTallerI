//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_JPL_EPH_DE430_H
#define PROYECTOTALLERI_JPL_EPH_DE430_H

#include "Matrix.h"
#include <vector>
using namespace std;

struct Ephemeris {
    Matrix r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter,
            r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon,
            r_Sun;

    Ephemeris()
            : r_Mercury(3,1),
              r_Venus  (3,1),
              r_Earth  (3,1),
              r_Mars   (3,1),
              r_Jupiter(3,1),
              r_Saturn (3,1),
              r_Uranus (3,1),
              r_Neptune(3,1),
              r_Pluto  (3,1),
              r_Moon   (3,1),
              r_Sun    (3,1)
    {}
};

Ephemeris JPL_Eph_DE430(double Mjd_TDB);



#endif //PROYECTOTALLERI_JPL_EPH_DE430_H
