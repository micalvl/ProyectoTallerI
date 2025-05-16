/**
 *  @file   G_AccelHarmonic.cpp
 *  @brief  G_AccelHarmonic function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-07
 ***********************************************/


/*
%--------------------------------------------------------------------------
%
% G_AccelHarmonic.m
%
% Purpose:
%   Computes the gradient of the Earth's harmonic gravity field
%
% Inputs:
%   r           Satellite position vector in the true-of-date system
%   U           Transformation matrix to body-fixed system
%   n           Gravity model degree
%   m 			Gravity model order
%
% Output:
%   G    		Gradient (G=da/dr) in the true-of-date system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/


#include "../include/G_AccelHarmonic.h"

Matrix G_AccelHarmonic(const Matrix& r, const Matrix& E, int n_max, int m_max) {
    const double d = 1.0;
    Matrix G(3, 3);
    Matrix dr(3, 1);
    for (int i = 1; i <= 3; ++i) {

        for (int k = 1; k <= 3; ++k)
            dr(k, 1) = 0.0;

        dr(i, 1) = d;

        Matrix r_plus  = r + dr.opsc(0.5);
        Matrix r_minus = r - dr.opsc(0.5);

        Matrix a_plus  = AccelHarmonic(r_plus, E, n_max, m_max);
        Matrix a_minus = AccelHarmonic(r_minus, E, n_max, m_max);

        Matrix da = a_plus - a_minus;

        for (int j = 1; j <= 3; ++j)
            G(j, i) = da(j, 1) / d;
    }

    return G;
}
