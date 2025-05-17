#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include "../include/Mjday.h"
#include "../include/Matrix.h"
#include "../include/R_x.h"
#include "../include/Legendre.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/gmst.h"
#include "../include/Frac.h"
#include "../include/AccelPointMass.h"
#include "../include/Sat_const.h"
#include "../include/Cheb3D.h"
#include "JPL_Eph_DE430.h"
#include "../include/MeasUpdate.h"
#include "MeanObliquity.h"
#include "AccelHarmonic.h"
#include "global.h"
#include "G_AccelHarmonic.h"
#include "sign_.h"
#include "Position.h"
#include "angl.h"
#include "NutAngles.h"
#include "EqnEquinox.h"

#define TOL_ 10e-14

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;

int Mjday_01()
{
    _assert(fabs(Mjday(2025,4,3,15,37,5)-60768.6507523148) < pow(10,-10));

    /*   cout << setprecision(20);
       cout << Mjday(2025,4,3,15,37,5) << endl;
       cout << 60768.6507523148;
      */
    return 0;
}

int Mjday_02()
{
    _assert(fabs(Mjday(2025,4,3,0,0,0.0)-Mjday(2025,4,3)) < pow(10,-10));

    return 0;
}

int R_x_01()
{
    double alpha = 1.0;
    Matrix sol(3, 3);

    sol = R_x(alpha);
    sol.print();

    _assert(fabs(sol(1,1)) - 1 < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3)) < TOL_);
    _assert(sol(2,1) < TOL_ && fabs(sol(2,2) - 0.54030230586814 ) < TOL_ && fabs(sol(2,3) -0.841470984807897) < TOL_);


    return 0;
}

int R_y_01()
{
    double alpha = 1.0;
    Matrix sol(3, 3);

    sol = R_y(alpha);
    sol.print();

    _assert(fabs(sol(1,1) - 0.54030230586814) < TOL_);
    _assert(fabs(sol(1,2)) < TOL_);
    _assert(fabs(sol(1,3) + 0.841470984807897) < TOL_);

    _assert(fabs(sol(2,1)) < TOL_);
    _assert(fabs(sol(2,2) - 1.0) < TOL_);
    _assert(fabs(sol(2,3)) < TOL_);

    _assert(fabs(sol(3,1) - 0.841470984807897) < TOL_);
    _assert(fabs(sol(3,2)) < TOL_);
    _assert(fabs(sol(3,3) - 0.54030230586814) < TOL_);

    return 0;
}

int R_z_01()
{
    double alpha = 1.0;
    Matrix sol(3, 3);

    sol = R_z(alpha);
    sol.print();

    _assert(fabs(sol(1,1) - 0.54030230586814) < TOL_);
    _assert(fabs(sol(1,2) - 0.841470984807897) < TOL_);
    _assert(fabs(sol(1,3)) < TOL_);

    _assert(fabs(sol(2,1) + 0.841470984807897) < TOL_);
    _assert(fabs(sol(2,2) - 0.54030230586814) < TOL_);
    _assert(fabs(sol(2,3)) < TOL_);

    _assert(fabs(sol(3,1)) < TOL_);
    _assert(fabs(sol(3,2)) < TOL_);
    _assert(fabs(sol(3,3) - 1.0) < TOL_);

    return 0;
}

int Legendre_01() {
    Matrix pnm(3, 3), dpnm(3, 3);
    double fi = 1.0;

    Legendre(2, 2, fi, pnm, dpnm);

    _assert(fabs(pnm(2, 2) - (sqrt(3.0) * cos(fi))) < TOL_);
    _assert(fabs(dpnm(2, 2) - (-sqrt(3.0) * sin(fi))) < TOL_);

    return 0;
}

int gmst_01() {
    double Mjd_24_04_2025 = 60825.0;
    double gmst_rad = gmst(Mjd_24_04_2025);
    double expected_gmst = 4.32424566600459492;

    _assert(fabs(gmst_rad - expected_gmst) < TOL_);
    return 0;
}

int frac_01(){
    _assert(fabs(Frac(1.5) - 0.5) < TOL_);
    _assert(fabs(Frac(2.25) - 0.25) < TOL_);
    _assert(fabs(Frac(1.66) - 0.66) < TOL_);
    return 0;
}

int meanObliquity01(){
    double t0 = MJD_J2000;
    double t1 = t0 + 36525.0;

    double m0 = MeanObliquity(t0);
    double m1 = MeanObliquity(t1);

    _assert(m1 < m0);

    return 0;
}

int meanObliquity02(){
    double Mjd_TT      = MJD_J2000 + 36525.0;
    double mobl_rad    = MeanObliquity(Mjd_TT);
    double expected    = 0.40886584462678882;

    _assert(fabs(mobl_rad - expected) < TOL_);

    return 0;
}

int measUpdate01(){
    Matrix x(1,1); x(1,1) = 0.0;
    Matrix z(1,1); z(1,1) = 4.0;
    Matrix g(1,1); g(1,1) = 2.0;
    Matrix s(1,1); s(1,1) = 2.0;
    Matrix G(1,1); G(1,1) = 1.0;
    Matrix P(1,1); P(1,1) = 1.0;

    MeasUpdate(x, z, g, s, G, P, 1);

    _assert(fabs(x(1,1) - 0.4) < TOL_);
    _assert(fabs(P(1,1) - 0.8) < TOL_);

    return 0;
}

int AccelHarmonic01() {

    Cnm = Matrix::zeros(301, 301);
    Snm = Matrix::zeros(201,  201);

    for (int n = 1; n < 300; ++n)
        for (int m = 1; m < 300; ++m)
            Cnm(n,m) = 0.0;
    for (int n = 1; n < 200; ++n)
        for (int m = 1; m < 200; ++m)
            Snm(n,m) = 0.0;

    Matrix r(3,1);
    r(1,1) = 7000e3;
    r(2,1) = 1234.5;
    r(3,1) = -2500.0;

    Matrix E = Matrix::identity(3);

    Matrix a = AccelHarmonic(r, E, 2, 2);

    _assert(fabs(a(1,1)) < TOL_);
    _assert(fabs(a(2,1)) < TOL_);
    _assert(fabs(a(3,1)) < TOL_);

    return 0;
}

int G_AccelHarmonic01() {
    Cnm = Matrix::zeros(301, 301);
    Snm = Matrix::zeros(201,  201);

    for (int n = 1; n < 300; ++n)
        for (int m = 1; m < 300; ++m)
            Cnm(n,m) = 0.0;
    for (int n = 1; n < 200; ++n)
        for (int m = 1; m < 200; ++m)
            Snm(n,m) = 0.0;
    Cnm(1,1) = 1.0;

    Matrix r(3,1);
    r(1,1) = 8000e3;
    r(2,1) = 1000e3;
    r(3,1) =  500e3;

    Matrix E = Matrix::identity(3);

    Matrix Gmat = G_AccelHarmonic(r, E, 0, 0);

    _assert(fabs(Gmat(1,2) - Gmat(2,1)) < TOL_);
    _assert(fabs(Gmat(1,3) - Gmat(3,1)) < TOL_);
    _assert(fabs(Gmat(2,3) - Gmat(3,2)) < TOL_);


    return 0;
}

int accel_point_mass_01()
{
    double r_vals[] = {2.0, 0.0, 0.0};
    double s_vals[] = {1.0, 0.0, 0.0};

    Matrix r(3, 1, r_vals, 3);
    Matrix s(3, 1, s_vals, 3);

    Matrix a = AccelPointMass(r, s, 1);

    _assert(fabs(a(1,1)+ 2.0) < TOL_);
    _assert(fabs(a(2,1)) < TOL_);
    _assert(fabs(a(3,1)) < TOL_);

    return 0;
}

int Cheb3D_01() {

    double Ta = 0.0, Tb = 1.0;
    int N = 1;
    double t = 0.5;

    vector<double> Cx = { 2.0 };
    vector<double> Cy = { 3.0 };
    vector<double> Cz = { 4.0 };

    Matrix result = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    _assert(fabs(result(1,1) - 2.0) < TOL_);
    _assert(fabs(result(1,2) - 3.0) < TOL_);
    _assert(fabs(result(1,3) - 4.0) < TOL_);

    return 0;
}

int sign_01() {
    _assert(fabs(sign_(-3.7,  1.0) - 3.7) < TOL_);
    _assert(fabs(sign_( 2.5, -1.0) + 2.5) < TOL_);
    return 0;
}
int position01() {
    Matrix r = Position(0.0, 0.0, 0.0);

    _assert(fabs(r(1,1) - R_Earth) < TOL_);
    _assert(fabs(r(2,1))           < TOL_);
    _assert(fabs(r(3,1))           < TOL_);

    return 0;
}

int angl01(){
    Matrix v1(2,1), v2(2,1);
    v1(1,1) = 1.0; v1(2,1) = 0.0;
    v2(1,1) = 0.0; v2(2,1) = 1.0;
    double angle1 = angl(v1, v2);

    _assert(fabs(angle1 - M_PI/2.0) < TOL_);

    return 0;
}

int angl02(){
    Matrix v1(2,1), v2(2,1);
    v1(1,1) = 1.0; v1(2,1) = 0.0;
    v2(1,1) = 1.0; v2(2,1) = 0.0;
    double angle = angl(v1, v2);

    _assert(fabs(angle - 0.0) < TOL_);

    return 0;
}

int angl03(){
    Matrix v1(2,1), v2(2,1);
    v1(1,1) = 0.0; v1(2,1) = 0.0;
    v2(1,1) = 1.0; v2(2,1) = 0.0;
    double angle = angl(v1, v2);

    _assert(fabs(angle - 999999.1) < TOL_);

    return 0;
}

/*
int JPL_Eph_DE430_01() {

    const char* path = "../data/DE430Coeff.txt";
    ifstream infile(path);
    _assert(infile.is_open());
    int rowCount = 0;
    string line;
    while (std::getline(infile, line)) {
        if (!line.empty()) ++rowCount;
    }
    infile.close();

    // 2) Definir el número de columnas que queremos leer
    const int colCount = 1020;

    // 3) Llenar la tabla PC
    DE430Coeff(rowCount, colCount);

    // 4) Elegir una fecha válida en el primer intervalo de PC
    double start = PC(1,1);
    double end   = PC(1,2);
    double Mjd_TDB = 0.5*(start + end);

    // 5) Llamar a JPL_Eph_DE430
    Ephemeris eph = JPL_Eph_DE430(Mjd_TDB);

    // 6) Comprobar dimensiones de salida
    _assert(eph.r_Earth.getFilas()   == 3);
    _assert(eph.r_Earth.getColumnas() == 1);
    _assert(eph.r_Sun.getFilas()     == 3);
    _assert(eph.r_Sun.getColumnas()   == 1);

    // 7) Comprobar normas razonables
    double normE = eph.r_Earth.norm();
    _assert(normE > 0.0 && normE < 5e11);

    double normS = eph.r_Sun.norm();
    _assert(normS > 1e10 && normS < 3e11);

    return 0;
}
*/

int nutangles01(){
    Matrix out = NutAngles(MJD_J2000);

    _assert(out.getFilas()   == 2);
    _assert(out.getColumnas()== 1);

    double dpsi = out(1,1);
    double deps = out(2,1);

    _assert(fabs(dpsi) < 0.01);
    _assert(fabs(deps) < 0.01);

}

int eqnEquinox01() {
    // 1) En J2000, la ecuación de los equinoccios es prácticamente cero
    double eq0 = EqnEquinox(MJD_J2000);
    _assert(fabs(eq0) < TOL_);

    return 0;
}



int all_tests()
{
    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(frac_01);
    _verify(gmst_01);
    _verify(Legendre_01);
    _verify(accel_point_mass_01);
    _verify(Cheb3D_01);
    //_verify(JPL_Eph_DE430_01);
    _verify(measUpdate01);
    _verify(meanObliquity01);
    _verify(meanObliquity02);
    _verify(AccelHarmonic01);
    _verify(G_AccelHarmonic01);
    _verify(sign_01);
    _verify(position01);
    _verify(angl01);
    _verify(angl02);
    _verify(angl03);
    _verify(nutangles01);
    _verify(eqnEquinox01);

    return 0;
}



int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);
/*
    double v[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {5.0, 6.0, 7.0, 8.0};

    Matrix m1(2, 2);
    m1.print();

    Matrix m2(2, 2, v, 4);
    m2.print();

    Matrix m3(2, 2, v2, 4);
    m3.print();

    m1 = m2 * m3 * m2;
    m1.print();
 */
    return result != 0;
}



