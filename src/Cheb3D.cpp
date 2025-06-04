/**
 *  @file   Cheb3D.cpp
 *  @brief  Cheb's function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-06
 ***********************************************/

#include "Cheb3D.h"


Matrix Cheb3D(double t, int N, double Ta, double Tb,
              const std::vector<double>& Cx,
              const std::vector<double>& Cy,
              const std::vector<double>& Cz)
{
    if (t < Ta || t > Tb)
        throw std::invalid_argument("Cheb3D: t out of [Ta,Tb]");

    double tau = 2.0*(t - Ta)/(Tb - Ta) - 1.0;

    Matrix f0(3,1), f1(3,1), f2(3,1);

    for (int k = N-1; k >= 1; --k) {
        f2 = f1;
        f1 = f0;

        Matrix coef(3,1);
        coef(1,1) = Cx[k];
        coef(2,1) = Cy[k];
        coef(3,1) = Cz[k];

        f0 = f1.opsc(2.0*tau) - f2 + coef;
    }

    Matrix C0(3,1);
    C0(1,1) = Cx[0];
    C0(2,1) = Cy[0];
    C0(3,1) = Cz[0];

    return f0.opsc(tau) - f1 + C0;
}

