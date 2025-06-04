/**
 *  @file   NutMatrix.cpp
 *  @brief  NutMatrix method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#include "../include/NutMatrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"

Matrix NutMatrix(double Mjd_TT) {
    double eps = MeanObliquity(Mjd_TT);
    Matrix nut = NutAngles(Mjd_TT);
    double dpsi = nut(1,1);
    double deps = nut(2,1);
    return R_x(-eps - deps) * R_z(-dpsi) * R_x(eps);
}
