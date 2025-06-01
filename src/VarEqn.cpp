//
// Created by micalvl on 03/04/2025.
//

#include "../include/VarEqn.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/AccelHarmonic.h"
#include "../include/G_AccelHarmonic.h"

Matrix VarEqn(double x, const Matrix& yPhi) {
    extern Param AuxParam;
    extern Matrix eopdata;

    IERSResult ier = IERS(eopdata, AuxParam.Mjd_UTC, 'l');
    double x_pole  = ier.x_pole;
    double y_pole  = ier.y_pole;
    double UT1_UTC = ier.UT1_UTC;
    double TAI_UTC = ier.TAI_UTC;
    double TT_UTC  = timediff(UT1_UTC, TAI_UTC).TT_UTC;
    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

    Matrix P = PrecMatrix(MJD_J2000, AuxParam.Mjd_TT + x/86400.0);
    Matrix N = NutMatrix(AuxParam.Mjd_TT + x/86400.0);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix r(3,1), v(3,1);
    for (int i=1; i<=3; ++i) {
        r(i,1) = yPhi(i,1);
        v(i,1) = yPhi(i+3,1);
    }

    Matrix Phi(6,6);
    for (int j=1; j<=6; ++j)
        for (int i=1; i<=6; ++i)
            Phi(i,j) = yPhi(6*(j-1)+i,1);

    Matrix a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
    Matrix G = G_AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

    Matrix dfdy = Matrix::zeros(6,6);
    for (int i=1; i<=3; ++i) {
        for (int j=1; j<=3; ++j) {
            dfdy(i,j)     = 0.0;
            dfdy(i+3,j)   = G(i,j);
            dfdy(i,j+3)   = (i==j) ? 1.0 : 0.0;
            dfdy(i+3,j+3) = 0.0;
        }
    }

    Matrix Phip = dfdy * Phi;

    Matrix yPhip(42,1);
    for (int i=1; i<=3; ++i) {
        yPhip(i,1) = v(i,1);
        yPhip(i+3,1) = a(i,1);
    }
    for (int j=1; j<=6; ++j)
        for (int i=1; i<=6; ++i)
            yPhip(6*(j-1)+i,1) = Phip(i,j);

    return yPhip;
}
