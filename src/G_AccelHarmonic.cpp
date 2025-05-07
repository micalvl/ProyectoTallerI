/**
 *  @file   G_AccelHarmonic.cpp
 *  @brief  G_AccelHarmonic function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-07
 ***********************************************/

#include "G_AccelHarmonic.h"
#include "Matrix.h"
#include "AccelHarmonic.h"

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


Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max) {
    double d = 1.0;
    Matrix G(3, 3);
    Matrix dr(3, 1);

    for (int i = 0; i < 3; ++i) {
        dr(0, 0) = 0.0;
        dr(1, 0) = 0.0;
        dr(2, 0) = 0.0;
        dr(i, 0) = d;


        Matrix r_plus  = r + numberMatrix(dr, 0.5);
        Matrix r_minus = r - (dr * 0.5);

        Matrix a_plus  = AccelHarmonic(r_plus, U, n_max, m_max);
        Matrix a_minus = AccelHarmonic(r_minus, U, n_max, m_max);
        Matrix da = a_plus - a_minus;

        for (int j = 0; j < 3; ++j) {
            G(j, i) = da(j, 0) / d;
        }
    }

    return G;
}
