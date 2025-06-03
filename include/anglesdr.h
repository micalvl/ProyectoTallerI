//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_ANGLESDR_H
#define PROYECTOTALLERI_ANGLESDR_H

#include "Matrix.h"

struct AnglesDRResult {
    Matrix r2;
    Matrix v2;
};

AnglesDRResult anglesdr(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        const Matrix& rsite1,
        const Matrix& rsite2,
        const Matrix& rsite3
);


#endif //PROYECTOTALLERI_ANGLESDR_H
