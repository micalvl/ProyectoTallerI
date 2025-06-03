//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_DEINTEG_H
#define PROYECTOTALLERI_DEINTEG_H


#include "Matrix.h"
#include <functional>
#include <vector>
#include <string>

enum class DE_STATE {
    DE_INIT = 1,
    DE_DONE = 2,
    DE_BADACC = 3,
    DE_NUMSTEPS = 4,
    DE_STIFF = 5,
    DE_INVPARAM = 6
};

using ODEFunction = std::function<Matrix(double, const Matrix&)>;

Matrix DEInteg(const ODEFunction& func, double t, double tout, double relerr,
               double abserr, int n_eqn, Matrix y);


#endif //PROYECTOTALLERI_DEINTEG_H
