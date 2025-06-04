/**
 *  @file   EccAnom.cpp
 *  @brief  EccAnom method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#include "../include/EqnEquinox.h"

double EqnEquinox(double Mjd_TT) {
    Matrix nut = NutAngles(Mjd_TT);
    double dpsi = nut(1,1);

    double eps = MeanObliquity(Mjd_TT);

    return dpsi * std::cos(eps);
}