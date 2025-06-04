/**
 *  @file   gibbs.cpp
 *  @brief  gibbs method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-14
 ***********************************************/

#ifndef PROYECTOTALLERI_GIBBS_H
#define PROYECTOTALLERI_GIBBS_H


#include "Matrix.h"
#include <string>
using namespace std;

struct GibbsResult {
    Matrix    v2;
    double    theta;
    double    theta1;
    double    copa;
    std::string error;
};

/// @param r1  3×1 position vector #1 [m]
/// @param r2  3×1 position vector #2 [m]
/// @param r3  3×1 position vector #3 [m]
/// @return    GibbsResult containing v2, angles and error flag
GibbsResult gibbs(const Matrix& r1,
                  const Matrix& r2,
                  const Matrix& r3);

#endif //PROYECTOTALLERI_GIBBS_H
