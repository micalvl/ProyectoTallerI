/**
 *  @file   GHAMatrix.cpp
 *  @brief  GHAMatrix method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#include "GHAMatrix.h"
#include "gast.h"
#include "R_z.h"
#include <cmath>

Matrix GHAMatrix(double Mjd_UT1) {
    double theta = gast(Mjd_UT1);
    return R_z(theta);
}
