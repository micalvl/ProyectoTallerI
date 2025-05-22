//
// Created by micalvl on 03/04/2025.
//

#include "PoleMatrix.h"
#include "R_x.h"
#include "R_y.h"

Matrix PoleMatrix(double xp, double yp) {
    return R_y(-xp) * R_x(-yp);
}
