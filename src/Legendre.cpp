/**
 *  @file   Legendre.cpp
 *  @brief  Legendre method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-10
 ***********************************************/

#include <cmath>
#include <algorithm>
#include "Legendre.h"



void Legendre(int n, int m, double phi, Matrix &pnm, Matrix &dpnm) {
    pnm = Matrix(n + 1, m + 1);
    dpnm = Matrix(n + 1, m + 1);

    double sinphi = sin(phi);
    double cosphi = cos(phi);

    pnm(1, 1) = 1.0;
    dpnm(1, 1) = 0.0;

    if (n == 0) return;

    pnm(2, 1) = sqrt(3.0) * sinphi;
    dpnm(2, 1) = sqrt(3.0) * cosphi;

    if (m >= 1) {
        pnm(2, 2) = sqrt(3.0) * cosphi;
        dpnm(2, 2) = -sqrt(3.0) * sinphi;
    }

    for (int i = 2; i <= n; ++i) {
        // Diagonal
        pnm(i + 1, i + 1) = sqrt((2.0 * i + 1.0) / (2.0 * i)) * cosphi * pnm(i, i);
        dpnm(i + 1, i + 1) = sqrt((2.0 * i + 1.0) / (2.0 * i)) *
                             (cosphi * dpnm(i, i) - sinphi * pnm(i, i));
    }

    for (int i = 1; i <= n; ++i) {
        pnm(i + 1, i) = sqrt(2.0 * i + 1.0) * sinphi * pnm(i, i);
        dpnm(i + 1, i) = sqrt(2.0 * i + 1.0) *
                         (cosphi * pnm(i, i) + sinphi * dpnm(i, i));
    }

    for (int j = 0; j <= m; ++j) {
        for (int i = j + 2; i <= n; ++i) {
            double a = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
            double b = sqrt(2.0 * i - 1.0);
            double c = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0));

            pnm(i + 1, j + 1) = a * (b * sinphi * pnm(i, j + 1) - c * pnm(i - 1, j + 1));
            dpnm(i + 1, j + 1) = a * (b * (sinphi * dpnm(i, j + 1) + cosphi * pnm(i, j + 1))
                                      - c * dpnm(i - 1, j + 1));
        }
    }
}


