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

Matrix Cheb3D(double t, int N, double Ta, double Tb,const vector<double>& Cx,const vector<double>& Cy,const vector<double>& Cz){
    //Check validity
    if( (t<Ta) || (Tb<t) ){
    cerr << "ERROR: Time out of range in Cheb3D::Value\n" << endl;
    }

    double tau = (2*t-Ta-Tb)/(Tb-Ta);

    Matrix f1(1, 3), f2(1, 3);

    for (int i = N; i >= 2; --i) {
        Matrix old_f1 = f1;

        double ci_vals[3] = {
                Cx[i-1],
                Cy[i-1],
                Cz[i-1]
        };
        Matrix Ci(1, 3, ci_vals, 3);

        f1 = f1.opsc(2.0*tau) - f2 + Ci;
        f2 = old_f1;
    }

    double c0_vals[3] = {
            Cx[0],
            Cy[0],
            Cz[0]
    };
    Matrix C0(1, 3, c0_vals, 3);

    Matrix ChebApp = f1.opsc(tau) - f2 + C0;

    return ChebApp;
}

