/**
 *  @file   Cheb3D.cpp
 *  @brief  Cheb's function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-06
 ***********************************************/

#include "Cheb3D.h"

/*
%--------------------------------------------------------------------------
%
% Chebyshev approximation of 3-dimensional vectors
%
% Inputs:
%     N       Number of coefficients
%     Ta      Begin interval
%     Tb      End interval
%     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
%     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
%     Cz      Coefficients of Chebyshev polyomial (z-coordinate)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
*/


Matrix Cheb3D(double t, int N, double Ta, double Tb,
              const std::vector<double>& Cx,
              const std::vector<double>& Cy,
              const std::vector<double>& Cz)
{
    if (t < Ta || t > Tb)
        throw std::invalid_argument("Cheb3D: t fuera de [Ta,Tb]");

    // tau en [-1,1]
    double tau = 2.0*(t - Ta)/(Tb - Ta) - 1.0;

    // f0, f1, f2 como vectores columna 3×1
    Matrix f0(3,1), f1(3,1), f2(3,1);

    // Clenshaw hacia atrás: k = N..1
    for (int k = N; k >= 1; --k) {
        // Desplazamos
        f2 = f1;
        f1 = f0;

        // Construimos coef(k) en 1-based
        Matrix coef(3,1);
        coef(1,1) = Cx[k];
        coef(2,1) = Cy[k];
        coef(3,1) = Cz[k];

        // f0 = 2·tau·f1 − f2 + coef
        f0 = f1.opsc(2.0*tau) - f2 + coef;
    }

    // Coeficiente de orden 0
    Matrix C0(3,1);
    C0(1,1) = Cx[0];
    C0(2,1) = Cy[0];
    C0(3,1) = Cz[0];

    // Evaluación final
    return f0.opsc(tau) - f1 + C0;
}

