//
// Created by micalvl on 03/04/2025.
//

#include "AzElPa.h"
#include <cmath>
using namespace std;

void AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds) {

    double s1 = s(1,1);
    double s2 = s(2,1);
    double s3 = s(3,1);

    double rho = sqrt(s1*s1 + s2*s2);

    Az = atan2(s1, s2);
    if (Az < 0.0) {
        Az += 2.0 * M_PI;
    }

    El = atan(s3 / rho);

    double rho2 = rho * rho;
    dAds = Matrix::zeros(3, 1);
    dAds(1,1) =  s2 / rho2;
    dAds(2,1) = -s1 / rho2;
    dAds(3,1) =  0.0;

    double denom = s1*s1 + s2*s2 + s3*s3;
    dEds = Matrix::zeros(3, 1);
    dEds(1,1) = (-s1 * s3 / rho) / denom;
    dEds(2,1) = (-s2 * s3 / rho) / denom;
    dEds(3,1) = ( rho           ) / denom;
}