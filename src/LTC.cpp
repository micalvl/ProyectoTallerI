//
// Created by micalvl on 03/04/2025.
//

#include "../include/LTC.h"
#include "R_y.h"
#include "R_z.h"

Matrix LTC(double lon, double lat) {
    Matrix M = R_y(-lat) * R_z(lon);
    for (int j = 1; j <= 3; ++j) {
        double Aux   = M(1,j);
        M(1,j) = M(2,j);
        M(2,j) = M(3,j);
        M(3,j) = Aux;
    }
    return M;
}