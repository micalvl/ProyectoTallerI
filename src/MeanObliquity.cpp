/**
 *  @file   MeanObliquity.cpp
 *  @brief  Computes the mean obliquity of the ecliptic from MJD_TT
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-07
 */


#include "Sat_const.h"
#include "../include/MeanObliquity.h"

/*
%--------------------------------------------------------------------------
%
% MeanObliquity.m
%
% Purpose:
%   Computes the mean obliquity of the ecliptic
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   MOblq     Mean obliquity of the ecliptic [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - MJD_J2000) / 36525.0;

    double arcsec = 84381.448 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T;

    return (3.14159265358979323846 / 180.0) * (arcsec / 3600.0);
}
