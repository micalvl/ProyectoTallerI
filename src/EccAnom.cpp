/**
 *  @file   EccAnom.cpp
 *  @brief  EccAnom method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-02
 ***********************************************/

#include "EccAnom.h"
#include <cmath>
#include <stdexcept>
#include <limits>
using namespace std;

double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 1;

    M = fmod(M, 2.0 * M_PI);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = M_PI;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    while (fabs(f) > 1e2 * numeric_limits<double>::epsilon()) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;
        if (i == maxit) {
            throw runtime_error("Convergence problems in EccAnom");
        }
    }

    return E;
}
