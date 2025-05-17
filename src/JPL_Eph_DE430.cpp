//
// Created by micalvl on 03/04/2025.
//

#include "JPL_Eph_DE430.h"
#include "Cheb3D.h"
#include "Sat_const.h"
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

/*
%--------------------------------------------------------------------------
%
% JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Inputs:
%   Mjd_TDB         Modified julian date of TDB
%
% Output:
%   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
%   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
%   r_Sun(geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Notes: Light-time is already taken into account
%
% Last modified:   2018/01/11   M. Mahooti
%
%--------------------------------------------------------------------------
*/

extern Matrix PC; // Problemas con esto

/// Construye coeficientes concatenando `blockCount` bloques de longitud
/// `blockLen`, arrancando en `startIndex` y saltando `step` casillas cada bloque.
static vector<double> getCoeffs(const vector<double>& pc,
                                     int startIndex,
                                     int blockCount,
                                     int step,
                                     int blockLen)
{
    vector<double> out;
    out.reserve(blockCount * blockLen);
    for (int b = 0; b < blockCount; ++b) {
        int base = startIndex - 1 + b*step;
        for (int k = 0; k < blockLen; ++k) {
            out.push_back(pc[base + k]);
        }
    }
    return out;
}

Ephemeris JPL_Eph_DE430(double Mjd_TDB) {
    Ephemeris eph;

    double JD = Mjd_TDB + 2400000.5;
    int rows = PC.getFilas(), cols = PC.getColumnas();
    int idx = 1;
    for (; idx <= rows; ++idx) {
        if (PC(idx,1) <= JD && JD <= PC(idx,2))
            break;
    }
    if (idx > rows)
        throw invalid_argument("JPL_Eph_DE430: Mjd_TDB fuera de rango en PC.");

    vector<double> PCrow(cols);
    for (int c = 1; c <= cols; ++c)
        PCrow[c-1] = PC(idx, c);

    double t1 = PCrow[0] - 2400000.5;
    double dt = Mjd_TDB - t1;

    // --- Earth ---
    {
        auto Cx = getCoeffs(PCrow, 231, 2, 39, 13);
        auto Cy = getCoeffs(PCrow, 244, 2, 39, 13);
        auto Cz = getCoeffs(PCrow, 257, 2, 39, 13);

        int blocks = 2;
        double interval = 16.0;
        int j = min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j*interval;

        Matrix v = Cheb3D(
                Mjd_TDB,       // t
                13,            // N
                Mjd0,          // Ta
                Mjd0 + interval,
                vector<double>(Cx.begin() + j*13, Cx.begin() + (j+1)*13),
                vector<double>(Cy.begin() + j*13, Cy.begin() + (j+1)*13),
                vector<double>(Cz.begin() + j*13, Cz.begin() + (j+1)*13)
        );
        eph.r_Earth = v.transpose().opsc(1e3);
    }

    // --- Moon ---
    {

        auto Cx = getCoeffs(PCrow, 441, 8, 39, 13);
        auto Cy = getCoeffs(PCrow, 454, 8, 39, 13);
        auto Cz = getCoeffs(PCrow, 467, 8, 39, 13);

        int blocks = 8;
        double interval = 4.0;
        int j = min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j*interval;

        Matrix v = Cheb3D(
                Mjd_TDB, 13, Mjd0, Mjd0 + interval,
                vector<double>(Cx.begin() + j*13, Cx.begin() + (j+1)*13),
                vector<double>(Cy.begin() + j*13, Cy.begin() + (j+1)*13),
                vector<double>(Cz.begin() + j*13, Cz.begin() + (j+1)*13)
        );
        eph.r_Moon = v.transpose().opsc(1e3);
    }

    // --- Sun ---
    {

        auto Cx = getCoeffs(PCrow, 753, 2, 33, 11);
        auto Cy = getCoeffs(PCrow, 764, 2, 33, 11);
        auto Cz = getCoeffs(PCrow, 775, 2, 33, 11);

        int blocks = 2;
        double interval = 16.0;
        int j = min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j*interval;

        Matrix v = Cheb3D(
                Mjd_TDB, 11, Mjd0, Mjd0 + interval,
                vector<double>(Cx.begin() + j*11, Cx.begin() + (j+1)*11),
                vector<double>(Cy.begin() + j*11, Cy.begin() + (j+1)*11),
                vector<double>(Cz.begin() + j*11, Cz.begin() + (j+1)*11)
        );
        eph.r_Sun = v.transpose().opsc(1e3);
    }

    // --- Mercury ---
    {

        auto Cx = getCoeffs(PCrow,   3, 4, 42, 14);
        auto Cy = getCoeffs(PCrow,  17, 4, 42, 14);
        auto Cz = getCoeffs(PCrow,  31, 4, 42, 14);

        int blocks = 4; double interval = 8.0;
        int j = min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j*interval;

        Matrix v = Cheb3D(
                Mjd_TDB, 14, Mjd0, Mjd0 + interval,
                vector<double>(Cx.begin() + j*14, Cx.begin() + (j+1)*14),
                vector<double>(Cy.begin() + j*14, Cy.begin() + (j+1)*14),
                vector<double>(Cz.begin() + j*14, Cz.begin() + (j+1)*14)
        );
        eph.r_Mercury = v.transpose().opsc(1e3);
    }

    // --- Venus ---
    {
        auto Cx = getCoeffs(PCrow, 171, 2, 30, 10);
        auto Cy = getCoeffs(PCrow, 181, 2, 30, 10);
        auto Cz = getCoeffs(PCrow, 191, 2, 30, 10);

        int blocks = 2; double interval = 16.0;
        int j = min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j*interval;

        Matrix v = Cheb3D(
                Mjd_TDB, 10, Mjd0, Mjd0 + interval,
                vector<double>(Cx.begin() + j*10, Cx.begin() + (j+1)*10),
                vector<double>(Cy.begin() + j*10, Cy.begin() + (j+1)*10),
                vector<double>(Cz.begin() + j*10, Cz.begin() + (j+1)*10)
        );
        eph.r_Venus = v.transpose().opsc(1e3);
    }

    // --- Mars ---
    {

        auto Cx = getCoeffs(PCrow, 309, 1,  0, 11);
        auto Cy = getCoeffs(PCrow, 320, 1,  0, 11);
        auto Cz = getCoeffs(PCrow, 331, 1,  0, 11);

        double interval = 32.0;
        Matrix v = Cheb3D(Mjd_TDB, 11, t1, t1 + interval, Cx, Cy, Cz);
        eph.r_Mars = v.transpose().opsc(1e3);
    }

    // --- Jupiter ---
    {
        auto Cx = getCoeffs(PCrow, 342, 1,  0,  8);
        auto Cy = getCoeffs(PCrow, 350, 1,  0,  8);
        auto Cz = getCoeffs(PCrow, 358, 1,  0,  8);

        double interval = 32.0;
        Matrix v = Cheb3D(Mjd_TDB, 8, t1, t1 + interval, Cx, Cy, Cz);
        eph.r_Jupiter = v.transpose().opsc(1e3);
    }

    // --- Saturn ---
    {
        auto Cx = getCoeffs(PCrow, 366, 1,  0, 7);
        auto Cy = getCoeffs(PCrow, 373, 1,  0, 7);
        auto Cz = getCoeffs(PCrow, 380, 1,  0, 7);

        double interval = 32.0;
        Matrix v = Cheb3D(Mjd_TDB, 7, t1, t1 + interval, Cx, Cy, Cz);
        eph.r_Saturn = v.transpose().opsc(1e3);
    }

    // --- Uranus ---
    {
        auto Cx = getCoeffs(PCrow, 387, 1,  0, 6);
        auto Cy = getCoeffs(PCrow, 393, 1,  0, 6);
        auto Cz = getCoeffs(PCrow, 399, 1,  0, 6);

        double interval = 32.0;
        Matrix v = Cheb3D(Mjd_TDB, 6, t1, t1 + interval, Cx, Cy, Cz);
        eph.r_Uranus = v.transpose().opsc(1e3);
    }

    // --- Neptune ---
    {
        auto Cx = getCoeffs(PCrow, 405, 1,  0, 6);
        auto Cy = getCoeffs(PCrow, 411, 1,  0, 6);
        auto Cz = getCoeffs(PCrow, 417, 1,  0, 6);

        double interval = 32.0;
        Matrix v = Cheb3D(Mjd_TDB, 6, t1, t1 + interval, Cx, Cy, Cz);
        eph.r_Neptune = v.transpose().opsc(1e3);
    }

    // --- Pluto ---
    {
        auto Cx = getCoeffs(PCrow, 423, 1,  0, 6);
        auto Cy = getCoeffs(PCrow, 429, 1,  0, 6);
        auto Cz = getCoeffs(PCrow, 435, 1,  0, 6);

        double interval = 32.0;
        Matrix v = Cheb3D(Mjd_TDB, 6, t1, t1 + interval, Cx, Cy, Cz);
        eph.r_Pluto = v.transpose().opsc(1e3);
    }


    double EMRAT  = 81.30056907419062;
    double EMRAT1 = 1.0 / (1.0 + EMRAT);
    eph.r_Earth   = eph.r_Earth - eph.r_Moon.opsc( EMRAT1 );
    eph.r_Mercury = eph.r_Mercury - eph.r_Earth;
    eph.r_Venus   = eph.r_Venus - eph.r_Earth;
    eph.r_Mars    = eph.r_Mars - eph.r_Earth;
    eph.r_Jupiter = eph.r_Jupiter - eph.r_Earth;
    eph.r_Saturn  = eph.r_Saturn - eph.r_Earth;
    eph.r_Uranus  = eph.r_Uranus - eph.r_Earth;
    eph.r_Neptune = eph.r_Neptune - eph.r_Earth;
    eph.r_Pluto   = eph.r_Pluto - eph.r_Earth;
    eph.r_Sun     = eph.r_Sun - eph.r_Earth;

    return eph;
}
