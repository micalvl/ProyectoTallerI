//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_ANGLESG_H
#define PROYECTOTALLERI_ANGLESG_H


#include "Matrix.h"

struct AnglesGResult {
    Matrix r2;
    Matrix v2;
};

AnglesGResult anglesg(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        const Matrix& Rs1, const Matrix& Rs2, const Matrix& Rs3
);

#endif //PROYECTOTALLERI_ANGLESG_H
