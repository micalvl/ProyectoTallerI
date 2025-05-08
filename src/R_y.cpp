/**
 *  @file   R_y.cpp
 *  @brief  R_y's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-09
 ***********************************************/

#include <valarray>
#include "../include/Matrix.h"
#include "../include/R_y.h"

/*
%--------------------------------------------------------------------------
%  input:
%    angle       - angle of rotation [rad]
%
%  output:
%    rotmat      - vector result
%--------------------------------------------------------------------------
 */

Matrix R_y(double alpha)
{
    Matrix rotmat(3,3);
    double C, S;

    C = cos(alpha);
    S = sin(alpha);
    // rotmat = zeros(3,3);

    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;


    return rotmat;
}
