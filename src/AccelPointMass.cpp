/**
 *  @file   AccelPointMass.cpp
 *  @brief  AccelPointMass's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-13
 ***********************************************/

#include "../include/AccelPointMass.h"


Matrix AccelPointMass(const Matrix& r, const Matrix& s, double GM) {

    if (r.getFilas() != 3 || r.getColumnas() != 1 ||
        s.getFilas() != 3 || s.getColumnas() != 1) {
        throw invalid_argument("r y s must be 3x1.");
    }

    Matrix d = r.operator-(s);

    double norm_d = d.norm();
    double norm_s = s.norm();


    Matrix term1 = d.opsc(1.0 / pow(norm_d, 3));
    Matrix term2 = s.opsc(1.0 / pow(norm_s, 3));

    Matrix a = (term1 + term2).opsc(-GM);

    return a;
}


