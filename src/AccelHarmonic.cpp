/**
 *  @file   AccelHarmonic.cpp
 *  @brief  AccelHarmonic function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-23
 ***********************************************/

#include "../include/AccelHarmonic.h"
#include "Legendre.h"
#include "Matrix.h"
#include <cmath>
#include <stdexcept>
#include "global.h"

using namespace std;

/*
%--------------------------------------------------------------------------
%
% AccelHarmonic.m
%
% Purpose:
%   Computes the acceleration due to the harmonic gravity field of the
%   central body
%
% Inputs:
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a           Acceleration (a=d^2r/dt^2)
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/

Matrix AccelHarmonic(const Matrix& r, const Matrix& E, int n_max, int m_max) {
    double r_ref = 6378.1363e3;
    double gm    = 398600.4415e9;

    if (Cnm.getFilas() < n_max+1 || Cnm.getColumnas() < m_max+1 ||
        Snm.getFilas() < n_max+1 || Snm.getColumnas() < m_max+1) {
        throw std::runtime_error("Cnm/Snm out of range");
    }

    Matrix r_bf = E * r;

    double d = r_bf.norm();
    double latgc = asin(r_bf(3,1) / d);
    double lon   = atan2(r_bf(2,1), r_bf(1,1));

    Matrix pnm(n_max+1, m_max+1);
    Matrix dpnm(n_max+1, m_max+1);

    Legendre(n_max, m_max, latgc, pnm, dpnm);

    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;

    for (int n=0; n<=n_max; ++n) {

        double q1 = 0.0, q2 = 0.0, q3 = 0.0;

        double b1 = (-gm / (d*d)) * pow(r_ref/d, n) * (n+1);
        double b2 =  (gm / d)     * pow(r_ref/d, n);
        double b3 =  (gm / d)     * pow(r_ref/d, n);

        for (int m=0; m<=m_max; ++m) {
            double coef_cnm = Cnm(n+1, m+1);
            double coef_snm = Snm(n+1, m+1);

            double cosmlon = cos(m * lon);
            double sinmlon = sin(m * lon);

            double common = coef_cnm * cosmlon + coef_snm * sinmlon;

            q1 += pnm(n+1, m+1) * common;
            q2 += dpnm(n+1, m+1) * common;
            q3 += m * pnm(n+1, m+1) * (coef_snm * cosmlon - coef_cnm * sinmlon);
        }

        dUdr     += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon   += q3 * b3;
    }

    double r2xy = pow(r_bf(1,1), 2) + pow(r_bf(2,1), 2);

    double ax = (1.0 / d * dUdr - r_bf(3,1) / (d*d*sqrt(r2xy)) * dUdlatgc) * r_bf(1,1)
                - (1.0 / r2xy * dUdlon) * r_bf(2,1);
    double ay = (1.0 / d * dUdr - r_bf(3,1) / (d*d*sqrt(r2xy)) * dUdlatgc) * r_bf(2,1)
                + (1.0 / r2xy * dUdlon) * r_bf(1,1);
    double az = (1.0 / d * dUdr) * r_bf(3,1) + (sqrt(r2xy) / (d*d)) * dUdlatgc;

    Matrix a_bf(3,1);
    a_bf(1,1) = ax;
    a_bf(2,1) = ay;
    a_bf(3,1) = az;

    Matrix a = E.transpose().operator*(a_bf);
    return a;
}



