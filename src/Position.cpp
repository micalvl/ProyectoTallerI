//
// Created by micalvl on 03/04/2025.
//

/*
%--------------------------------------------------------------------------
%
% Position.m
%
% Purpose:
%   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
%   latitude [rad], altitude [m])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
#include "../include/Position.h"



Matrix Position(double lon, double lat, double h) {

    double e2 = f_Earth * (2.0 - f_Earth);

    double CosLat = cos(lat);
    double SinLat = sin(lat);

    double N = R_Earth / sqrt(1.0 - e2 * SinLat * SinLat);

    Matrix r(3, 1);
    r(1,1) = (N + h) * CosLat * cos(lon);
    r(2,1) = (N + h) * CosLat * sin(lon);
    r(3,1) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}