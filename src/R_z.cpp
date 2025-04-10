//
// Created by micalvl on 03/04/2025.
//

#include <cmath>
#include "../include/R_z.h"

/*
%--------------------------------------------------------------------------
%  input:
%    angle       - angle of rotation [rad]
%
%  output:
%    rotmat      - vector result
%--------------------------------------------------------------------------
*/


Matrix R_z(double alpha) {
    Matrix rotmat(3, 3);
    double C, S;

    C = cos(alpha);
    S = sin(alpha);
    //rotmat = zeros(3,3);

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