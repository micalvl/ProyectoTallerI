//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_HGIBBS_H
#define PROYECTOTALLERI_HGIBBS_H


#include "Matrix.h"
#include <string>
#include "Sat_const.h"
#include "angl.h"
#include "unit.h"

struct HGibbsResult {
    Matrix    v2;
    double    theta;
    double    theta1;
    double    copa;
    std::string error;
};

HGibbsResult hgibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3, double Mjd1, double Mjd2, double Mjd3);

#endif //PROYECTOTALLERI_HGIBBS_H
