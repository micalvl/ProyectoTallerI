/**
 *  @file   Cheb3D.h
 *  @brief  Cheb's function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-06
 ***********************************************/

#ifndef PROYECTOTALLERI_CHEB3D_H
#define PROYECTOTALLERI_CHEB3D_H

#include <iostream>
#include "../include/Matrix.h"
#include <vector>
using namespace std;


/**
 * @brief Chebyshev approximation of 3-dimensional vectors.
 * @param[in] t Time.
 * @param[in] N Number of coefficients
 * @param[in] Ta Begin interval.
 * @param[in] Tb End interval.
 * @param[in] Cx Coefficients of Chebyshev polyomial (x-coordinate).
 * @param[in] Cy Coefficients of Chebyshev polyomial (y-coordinate).
 * @param[in] Cz Coefficients of Chebyshev polyomial (z-coordinate).
 * @return Matrix 3×1 column vector containing the interpolated (x, y, z) values at time t.
 */
Matrix Cheb3D(double t, int N, double Ta, double Tb, const vector<double>& Cx, const vector<double>& Cy, const vector<double>& Cz);


#endif //PROYECTOTALLERI_CHEB3D_H
