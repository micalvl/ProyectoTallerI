/**
 *  @file   PoleMatrix.cpp
 *  @brief  PoleMatrix method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#include "PoleMatrix.h"
#include "R_x.h"
#include "R_y.h"

Matrix PoleMatrix(double xp, double yp) {
    return R_y(-xp) * R_x(-yp);
}
