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

vector<double> Cheb3D(double t, int N, double Ta, double Tb,const vector<double>& Cx,const vector<double>& Cy,const vector<double>& Cz){
    //Check validity
    if( (t<Ta) || (Tb<t) ){
    cerr << "ERROR: Time out of range in Cheb3D::Value\n" << endl;
    }

    double tau = (2*t-Ta-Tb)/(Tb-Ta);

    vector<double> f1(3, 0.0);
    vector<double> f2(3, 0.0);

    for (int i = N - 1; i >= 1; --i) {

        vector<double> old_f1 = f1;

        f1[0] = 2 * tau * f1[0] - f2[0] + Cx[i];
        f1[1] = 2 * tau * f1[1] - f2[1] + Cy[i];
        f1[2] = 2 * tau * f1[2] - f2[2] + Cz[i];

        f2 = old_f1;
    }

    vector<double> ChebApp(3, 0.0);
    ChebApp[0] = tau * f1[0] - f2[0] + Cx[0];
    ChebApp[1] = tau * f1[1] - f2[1] + Cy[0];
    ChebApp[2] = tau * f1[2] - f2[2] + Cz[0];

    return ChebApp;
}

