/**
 *  @file   DEInteg.h
 *  @brief  DeInteg method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

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


/**
 * @brief Numerical integration methods for ordinaray differential equations. This module provides implemenation
 * of the variable order variable stepsize multistep method of Shampine & Gordon.
 * @param[in] func ODEFunction function.
 * @param[in] t Initial time value.
 * @param[in] tout Final time value.
 * @param[in] relerr Relative error tolerance.
 * @param[in] abserr Absolute error tolerance.
 * @param[in] n_eqn Number of equations.
 * @param[in,out] y DEInteg matrix.
 * @return Matrix y modified.
 */
Matrix DEInteg(const ODEFunction& func, double t, double tout, double relerr,
               double abserr, int n_eqn, Matrix y);


#endif //PROYECTOTALLERI_DEINTEG_H
