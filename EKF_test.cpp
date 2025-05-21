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
#include "gibbs.h"
#include "unit.h"
#include "hgibbs.h"
#include "gast.h"
#include "elements.h"
#include "IERS.h"
#include "LTC.h"
#include "timediff.h"
#include "Geodetic.h"

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



int JPL_Eph_DE430_01() {
    // 1) Verificar fichero
    const char* path = "../data/DE430Coeff.txt";
    std::ifstream f(path);
    _assert(f.is_open());
    f.close();

    // 2) Cargar PC en JD
    const int ROWS = 2285, COLS = 1020;
    DE430Coeff(ROWS, COLS);

    // 3) Imprimir primer intervalo (JD)
    double JD0 = PC(1,1), JD1 = PC(1,2);
    std::cout << "PC(1,1) = " << JD0 << ", PC(1,2) = " << JD1 << "  (JD)\n";

    // 4) Calcular MJD_TDB
    double Mjd_TDB = 0.5 * (JD0 + JD1) - 2400000.5;
    _assert(Mjd_TDB > 0.0);
    std::cout << "Mjd_TDB = " << Mjd_TDB << "\n";

    // 5) Llamar a jugadora principal
    Ephemeris eph = JPL_Eph_DE430(Mjd_TDB);

    // 6) Corregir orientación si es 1×3
    if (eph.r_Earth.getFilas()==1 && eph.r_Earth.getColumnas()==3)
        eph.r_Earth = eph.r_Earth.transpose();
    if (eph.r_Sun.getFilas()==1 && eph.r_Sun.getColumnas()==3)
        eph.r_Sun   = eph.r_Sun.transpose();

    // 7) Verificar dimensiones
    _assert(eph.r_Earth.getFilas()==3 && eph.r_Earth.getColumnas()==1);
    _assert(eph.r_Sun.getFilas()==3   && eph.r_Sun.getColumnas()==1);

    // 8) Verificar normas razonables
    double ne = eph.r_Earth.norm();
    double ns = eph.r_Sun.norm();
    _assert(ne > TOL_);
    _assert(ns > TOL_);

    std::cout << "[OK] Test JPL_Eph_DE430_01 superado.\n";
    return 0;
}


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
    double eq0 = EqnEquinox(MJD_J2000);

    _assert(fabs(eq0) < 0.01);

    return 0;
}

int gibbs01() {
    using namespace std;

    Matrix r1(3,1), r2(3,1), r3(3,1);
    r1(1,1)=1;  r1(2,1)=0;  r1(3,1)=0;
    r2(1,1)=0;  r2(2,1)=1;  r2(3,1)=0;
    r3(1,1)=-1; r3(2,1)=0;  r3(3,1)=0;

    GibbsResult g = gibbs(r1, r2, r3);
    _assert(g.error == "ok");
    _assert(fabs(g.theta  - M_PI/2)  < TOL_);
    _assert(fabs(g.theta1 - M_PI/2)  < TOL_);
    _assert(fabs(g.copa)             < TOL_);
    _assert(g.v2.norm() > 0.0);

    return 0;
}

int unit01() {
    Matrix v(3,1);
    v(1,1) = 1.0; v(2,1) = 0.0; v(3,1) = 0.0;
    Matrix u = unit(v);
    _assert(u(1,1) == 1.0);
    _assert(u(2,1) == 0.0);
    _assert(u(3,1) == 0.0);
    return 0;
}

int hgibbs01() {
    Matrix r1(3,1), r2(3,1), r3(3,1);
    double a1 = 0.01;   // ≈0.57°
    double a2 = 0.015;  // ≈0.86°
    r1(1,1)=1; r1(2,1)=0; r1(3,1)=0;
    r2(1,1)=std::cos(a1); r2(2,1)=std::sin(a1); r2(3,1)=0;
    r3(1,1)=std::cos(a2); r3(2,1)=std::sin(a2); r3(3,1)=0;

    HGibbsResult res = hgibbs(r1, r2, r3, 0.0, 1.0, 2.0);
    const double tol_angle = 0.01745329251994; // The same as the function

    _assert(res.error == "ok");
    _assert(res.theta  > 0.0 && res.theta  < tol_angle);
    _assert(res.theta1 > 0.0 && res.theta1 < tol_angle);
    _assert(fabs(res.copa) < TOL_);
    _assert(res.v2.norm() > 0.0);

    return 0;
}

int gast01() {
    double Mjd = MJD_J2000;
    double g1  = gast(Mjd);
    double g2  = gmst(Mjd) + EqnEquinox(Mjd);
    g2 = std::fmod(g2, 2*M_PI);
    if (g2 < 0) g2 += 2*M_PI;

    _assert(fabs(g1 - g2) < TOL_);

    return 0;
}


int elements01() {
    double vals[6] = {1.0, 1.0, 1.0,
                      1.0, 1.0, 2.0};

    Matrix y(6,1, vals, 6);
    ElementsResult el = elements(y);

    _assert(fabs(el.p - el.a*(1.0 - el.e*el.e)) < TOL_);

    return 0;
}

int IERS01() {
    double Mjd_UTC = 30000.0;

    Matrix eop(13,1);
    for (int r = 1; r <= 13; ++r) {
        eop(r,1) = 0.0;
    }
    eop(4,1) = Mjd_UTC;

    //  3600 = π/180 rad
    eop(5,1)  =  3600;
    eop(6,1)  =  7200;
    eop(7,1)  =   3.5;
    eop(8,1)  =   0.1;
    eop(9,1)  = 10800;
    eop(10,1) = 14400;
    eop(11,1) =  1800;
    eop(12,1) =  2700;
    eop(13,1) =   37.0;

    IERSResult R = IERS(eop, Mjd_UTC, 'n');

    _assert(fabs(R.x_pole   - (M_PI/180.0))        < TOL_);
    _assert(fabs(R.y_pole   - (2.0*M_PI/180.0))    < TOL_);
    _assert(fabs(R.UT1_UTC  - 3.5)                 < TOL_);
    _assert(fabs(R.LOD      - 0.1)                 < TOL_);
    _assert(fabs(R.dpsi     - (3.0*M_PI/180.0))    < TOL_);
    _assert(fabs(R.deps     - (4.0*M_PI/180.0))    < TOL_);
    _assert(fabs(R.dx_pole  - (0.5*M_PI/180.0))    < TOL_);
    _assert(fabs(R.dy_pole  - (0.75*M_PI/180.0))   < TOL_);
    _assert(fabs(R.TAI_UTC  - 37.0)                < TOL_);

    return 0;
}

int LTC01() {
    Matrix M = LTC(0.0, 0.0);

    _assert(fabs(M(1,1) - 0.0) < TOL_);
    _assert(fabs(M(1,2) - 1.0) < TOL_);
    _assert(fabs(M(1,3) - 0.0) < TOL_);
    _assert(fabs(M(2,1) - 0.0) < TOL_);
    _assert(fabs(M(2,2) - 0.0) < TOL_);
    _assert(fabs(M(2,3) - 1.0) < TOL_);
    _assert(fabs(M(3,1) - 1.0) < TOL_);
    _assert(fabs(M(3,2) - 0.0) < TOL_);
    _assert(fabs(M(3,3) - 0.0) < TOL_);

    return 0;
}

int timediff01(){
    TimeDiffResult R = timediff(1.0, 2.0);

    _assert(fabs(R.UT1_TAI - (-1.0))   < TOL_);
    _assert(fabs(R.UTC_GPS - 17.0)     < TOL_);
    _assert(fabs(R.UT1_GPS - 18.0)     < TOL_);
    _assert(fabs(R.TT_UTC  - 34.184)   < TOL_);
    _assert(fabs(R.GPS_UTC - (-17.0))  < TOL_);

    return 0;
}

int Geodetic01() {

    Matrix r(3,1);
    r(1,1) = R_Earth + 1000.0;
    r(2,1) = 0.0;
    r(3,1) = 0.0;

    double lon, lat, h;
    Geodetic(r, lon, lat, h);

    _assert(fabs(lon) < TOL_);
    _assert(fabs(lat) < TOL_);
    _assert(fabs(h - 1000.0) < TOL_);

    return 0;
}




    int all_tests()
{
    _verify(JPL_Eph_DE430_01);
    /*
    _verify(AccelHarmonic01);
    _verify(G_AccelHarmonic01);
    _verify(sign_01);
    _verify(position01);
    _verify(angl01);
    _verify(angl02);
    _verify(angl03);
    //_verify(nutangles01);
    _verify(eqnEquinox01);
    _verify(gast01);
    _verify(elements01);
    _verify(IERS01);
    _verify(LTC01);
    _verify(timediff01);
    _verify(gibbs01);
    _verify(unit01);
    _verify(hgibbs01);
    _verify(Legendre_01);
    _verify(accel_point_mass_01);
    _verify(Cheb3D_01);


    _verify(measUpdate01);
    _verify(meanObliquity01);
    _verify(meanObliquity02);


    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(frac_01);
    _verify(gmst_01);
    _verify(Geodetic01);
*/

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



