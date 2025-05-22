//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_DOUBLER_H
#define PROYECTOTALLERI_DOUBLER_H

#include "Matrix.h"
#include "Sat_const.h"


struct DoublerResult {
    Matrix r2, r3;
    double f1, f2, q1;
    double magr1, magr2, a, deltae32;
};

DoublerResult doubler(
        double cc1, double cc2,
        double magrsite1, double magrsite2,
        double magr1in,   double magr2in,
        const Matrix& los1, const Matrix& los2, const Matrix& los3,
        const Matrix& rsite1,const Matrix& rsite2,const Matrix& rsite3,
        double t1, double t3,
        bool direct
);


#endif //PROYECTOTALLERI_DOUBLER_H
