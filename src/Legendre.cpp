/**
 *  @file   Legendre.cpp
 *  @brief  Legendre method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-10
 ***********************************************/

#include <cmath>
#include "Legendre.h"



void Legendre(int n, int m, double fi, Matrix &pnm, Matrix &dpnm) {
    pnm = Matrix(n + 1, m + 1);
    dpnm = Matrix(n + 1, m + 1);

    pnm(1, 1) = 1.0;
    dpnm(1, 1) = 0.0;

    if (n >= 1 && m >= 1) {
        double sq3 = sqrt(3.0);
        pnm(2,2)  = sq3 * cos(fi);
        dpnm(2,2) = -sq3 * sin(fi);
    }

    for (int i = 2; i <= n; i++) {
        pnm(i + 1, i + 1) = sqrt((2.0 * i + 1.0) / (2.0 * i)) * cos(fi) * pnm(i, i);
        dpnm(i + 1, i + 1) = sqrt((2.0 * i + 1.0) / (2.0 * i)) * (cos(fi) * dpnm(i, i) - sin(fi) * pnm(i, i));
    }

    for (int i = 1; i <= n; i++) {
        pnm(i + 1, i) = sqrt(2.0 * i + 1.0) * sin(fi) * pnm(i, i);
        dpnm(i + 1, i) = sqrt(2.0 * i + 1.0) * (cos(fi) * pnm(i, i) + sin(fi) * dpnm(i, i));
    }

    for (int j = 0; j <= m; j++) {
        for (int i = j + 2; i <= n; i++) {
            double denom = (i - j) * (i + j);
            if (denom == 0.0) continue;

            double a = sqrt((2.0 * i + 1.0) / denom);
            double b = sqrt(2.0 * i - 1.0);
            double c = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0));

            pnm(i + 1, j + 1) = a * (b * sin(fi) * pnm(i, j + 1) - c * pnm(i - 1, j + 1));
            dpnm(i + 1, j + 1) = a * (b * (sin(fi) * dpnm(i, j + 1) + cos(fi) * pnm(i, j + 1)) - c * dpnm(i - 1, j + 1));
        }
    }
}

