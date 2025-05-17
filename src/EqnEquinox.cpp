//
// Created by micalvl on 03/04/2025.
//

#include "EqnEquinox.h"

double EqnEquinox(double Mjd_TT) {
    Matrix nut = NutAngles(Mjd_TT);
    double dpsi = nut(1,1);

    double eps = MeanObliquity(Mjd_TT);

    return dpsi * std::cos(eps);
}