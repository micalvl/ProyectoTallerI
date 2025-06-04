/**
 *  @file   R_z.cpp
 *  @brief  R_z's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-09
 ***********************************************/

#include <cmath>
#include "../include/R_z.h"




Matrix R_z(double angle) {
    double C, S;

    C = cos(angle);
    S = sin(angle);
    Matrix rotmat = Matrix::zeros(3,3);

    rotmat(1, 1) = C;
    rotmat(1, 2) = S;
    rotmat(1, 3) = 0.0;
    rotmat(2, 1) = -1.0 * S;
    rotmat(2, 2) = C;
    rotmat(2, 3) = 0.0;
    rotmat(3, 1) = 0.0;
    rotmat(3, 2) = 0.0;
    rotmat(3, 3) = 1.0;

    return rotmat;
}