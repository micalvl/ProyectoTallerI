/**
 *  @file   AccelPointMass.cpp
 *  @brief  AccelPointMass's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-13
 ***********************************************/

#include "../include/AccelPointMass.h"

/*
%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%				  mass
%
% Inputs:
%   r           Satellite position vector
%   s           Point mass position vector
%   GM          Gravitational coefficient of point mass
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%---------------------------------------------------------------------
 */



array<double, 3> AccelPointMass(const array<double, 3> r, const array<double, 3> s, double GM) {

    array<double, 3> a;

    double d0 = r[0] - s[0];
    double d1 = r[1] - s[1];
    double d2 = r[2] - s[2];

    double norm_d = sqrt(d0*d0 + d1*d1 + d2*d2);
    double norm_s = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);

    double norm_d3 = pow(norm_d, 3);
    double norm_s3 = pow(norm_s, 3);

    // a = -GM * ( d/norm(d)^3 + s/norm(s)^3 )
    a[0] = -GM * (d0 / norm_d3 + s[0] / norm_s3);
    a[1] = -GM * (d1 / norm_d3 + s[1] / norm_s3);
    a[2] = -GM * (d2 / norm_d3 + s[2] / norm_s3);

    return a;
}


