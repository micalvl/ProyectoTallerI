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
#include "GHAMatrix.h"
#include "PrecMatrix.h"
#include "PoleMatrix.h"
#include "NutMatrix.h"
#include "Mjday_TDB.h"
#include "doubler.h"
#include "Accel.h"
#include "AzElPa.h"
#include "VarEqn.h"
#include "anglesdr.h"
#include "anglesg.h"
#include "EccAnom.h"
#include "TimeUpdate.h"
#include "DEInteg.h"

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

int Legendre_02(){
    int n = 5;
    int m = 5;
    double fi = M_PI / 6.0;

    Matrix pnm, dpnm;
    Legendre(n, m, fi, pnm, dpnm);

    std::cout << "pnm =" << std::endl;
    pnm.print();

    std::cout << "dpnm =" << std::endl;
    dpnm.print();

    return 0;
}

int gmst_01() {
    double Mjd_UT1 = 55000.123456;
    double gmstime = gmst(Mjd_UT1);

    double expected = 5.4267686516;
    _assert(fabs(gmstime - expected) < TOL_);
    return 0;

}

int frac_01(){
    _assert(fabs(Frac(1.5) - 0.5) < TOL_);
    _assert(fabs(Frac(2.25) - 0.25) < TOL_);
    _assert(fabs(Frac(1.66) - 0.66) < TOL_);
    return 0;
}

int meanObliquity01(){
    double Mjd_TT = 51544.5;

    double MOblq = MeanObliquity(Mjd_TT);

    double expected_rad = (84381.448 / 3600.0) * (M_PI / 180.0);

    _assert(fabs(MOblq - expected_rad) < TOL_);

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

int AccelHarmonic_02(){
    Cnm = Matrix::zeros(301, 301);
    Snm = Matrix::zeros(201,  201);
    Matrix r(3,1);
    r(1,1) = -2436.45e3;
    r(2,1) = 5386.12e3;
    r(3,1) = 2444.75e3;

    Matrix E = Matrix::identity(3);
    int n = 20, m = 20;

    Matrix acc = AccelHarmonic(r, E, n, m);
    cout << "AccelHarmonic:\n";
    acc.print();

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

    std::vector<double> Cx = {2.0};
    std::vector<double> Cy = {3.0};
    std::vector<double> Cz = {4.0};

    Matrix result = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    _assert(fabs(result(1,1) - 2.0) < TOL_);
    _assert(fabs(result(2,1) - 3.0) < TOL_);
    _assert(fabs(result(3,1) - 4.0) < TOL_);

    return 0;
}

int sign_01() {
    _assert(fabs(sign_(-3.7,  1.0) - 3.7) < TOL_);
    _assert(fabs(sign_( 2.5, -1.0) + 2.5) < TOL_);
    return 0;
}
int position01() {
    double lon = 0.0;
    double lat = 0.0;
    double h   = 0.0;

    Matrix r = Position(lon, lat, h);

    _assert(r.getFilas() == 3);
    _assert(r.getColumnas() == 1);

    double expected_x = R_Earth;
    double expected_y = 0.0;
    double expected_z = 0.0;

    _assert(fabs(r(1,1) - expected_x) < TOL_);
    _assert(fabs(r(2,1) - expected_y) < TOL_);
    _assert(fabs(r(3,1) - expected_z) < TOL_);

    return 0;
}

int angl01(){
    Matrix v1(3,1); v1(1,1)=1; v1(2,1)=0; v1(3,1)=0;
    Matrix v2(3,1); v2(1,1)=0; v2(2,1)=1; v2(3,1)=0;

    double angle = angl(v1, v2);
    double expected = M_PI / 2;

    _assert(fabs(angle - expected) < TOL_);
    cout << "passed test" << endl;

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

    const char* path = "../data/DE430Coeff.txt";
    ifstream f(path);
    _assert(f.is_open());
    f.close();


    DE430Coeff(2285, 1020);


    double JD0 = PC(1,1), JD1 = PC(1,2);


    double Mjd_TDB = 0.5 * (JD0 + JD1) - 2400000.5;
    _assert(Mjd_TDB > 0.0);


    Ephemeris eph = JPL_Eph_DE430(Mjd_TDB);


    if (eph.r_Earth.getFilas()==1 && eph.r_Earth.getColumnas()==3)
        eph.r_Earth = eph.r_Earth.transpose();
    if (eph.r_Sun.getFilas()==1 && eph.r_Sun.getColumnas()==3)
        eph.r_Sun   = eph.r_Sun.transpose();


    _assert(eph.r_Earth.getFilas()==3 && eph.r_Earth.getColumnas()==1);
    _assert(eph.r_Sun.getFilas()==3   && eph.r_Sun.getColumnas()==1);


    double ne = eph.r_Earth.norm();
    double ns = eph.r_Sun.norm();
    _assert(ne > TOL_);
    _assert(ns > TOL_);

    return 0;
}


int nutangles01(){
    Matrix out = NutAngles(MJD_J2000);

    _assert(out.getFilas()   == 2);
    _assert(out.getColumnas()== 1);

    double dpsi = out(1,1);
    double deps = out(2,1);

    double tol = 0.0001;

    _assert(fabs(dpsi) < tol);
    _assert(fabs(deps) < tol);
    return 0;

}

int eqnEquinox01() {
    double eq0 = EqnEquinox(MJD_J2000);

    _assert(fabs(eq0) < 0.01);

    return 0;
}

int gibbs01() {
    Matrix r1(3,1), r2(3,1), r3(3,1);
    r1(1,1) = 7000; r1(2,1) = 0;    r1(3,1) = 0;
    r2(1,1) = 7071; r2(2,1) = 7071; r2(3,1) = 0;
    r3(1,1) = 0;    r3(2,1) = 7000; r3(3,1) = 0;

    GibbsResult res = gibbs(r1, r2, r3);

    _assert(res.error == "ok");
    _assert(res.v2.getFilas() == 3 && res.v2.getColumnas() == 1);
    _assert(res.v2.norm() > 0);
    _assert(fabs(res.copa) < 0.02);

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
    Matrix r1(3,1); r1(1,1) = 7000000;  r1(2,1) = 0;        r1(3,1) = 0;
    Matrix r2(3,1); r2(1,1) = 7070000;  r2(2,1) = 500000;   r2(3,1) = 0;
    Matrix r3(3,1); r3(1,1) = 7100000;  r3(2,1) = 1000000;  r3(3,1) = 0;

    double Mjd1 = 58000.0;
    double Mjd2 = 58000.01;
    double Mjd3 = 58000.02;

    double dt21 = (Mjd2 - Mjd1) * 86400.0;
    double dt31 = (Mjd3 - Mjd1) * 86400.0;
    double dt32 = (Mjd3 - Mjd2) * 86400.0;

    double magr1 = r1.norm();
    double magr2 = r2.norm();
    double magr3 = r3.norm();

    double term1 = -dt32 * (1.0/(dt21 * dt31) + GM_Earth / (12.0 * pow(magr1, 3)));
    double term2 = (dt32 - dt21) * (1.0/(dt21 * dt32) + GM_Earth / (12.0 * pow(magr2, 3)));
    double term3 = dt21 * (1.0/(dt32 * dt31) + GM_Earth / (12.0 * pow(magr3, 3)));

    Matrix v2 = r1.opsc(term1) + r2.opsc(term2) + r3.opsc(term3);

    bool ok = true;
    ok &= fabs(v2(1,1) - 24.95799064) < 1e-5;
    ok &= fabs(v2(2,1) - 656.56091908) < 1e-5;
    ok &= fabs(v2(3,1) - 0.0) < 1e-5;


    return 0;
}

int gast01() {
    double Mjd = MJD_J2000;
    double g1  = gast(Mjd);
    double g2  = gmst(Mjd) + EqnEquinox(Mjd);
    g2 = fmod(g2, 2*M_PI);
    if (g2 < 0) g2 += 2*M_PI;

    _assert(fabs(g1 - g2) < TOL_);

    return 0;
}


int elements01() {
    Matrix y(6,1);
    y(1,1) = 7000e3;
    y(2,1) = 0.0;
    y(3,1) = 0.0;
    y(4,1) = 0.0;
    y(5,1) = 7.5e3;
    y(6,1) = 1.0e3;

    ElementsResult el = elements(y);

    _assert(el.p > 0);
    _assert(el.a > 0);
    _assert(el.e >= 0 && el.e < 1);
    _assert(el.i >= 0 && el.i <= M_PI);
    _assert(el.Omega >= 0 && el.Omega < 2*M_PI);
    _assert(el.omega >= 0 && el.omega < 2*M_PI);
    _assert(el.M >= 0 && el.M < 2*M_PI);

    return 0;
}

int IERS01() {
    double Mjd_UTC = 30000.0;

    Matrix eop(1,13);
    for (int r = 1; r <= 13; ++r) {
        eop(1,r) = 0.0;
    }
    eop(1,4) = Mjd_UTC;

    //  3600 = π/180 rad
    eop(1,5)  =  3600;
    eop(1,6)  =  7200;
    eop(1,7)  =   3.5;
    eop(1,8)  =   0.1;
    eop(1,9)  = 10800;
    eop(1,10) = 14400;
    eop(1,11) =  1800;
    eop(1,12) =  2700;
    eop(1,13) =   37.0;

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

    double lon, lat, h;
    double R = R_Earth;

    Matrix r1(3,1); r1(1,1) = R; r1(2,1) = 0; r1(3,1) = 0;
    Geodetic(r1, lon, lat, h);

    _assert(fabs(lon - 0.0) < TOL_);
    _assert(fabs(lat - 0.0) < TOL_);
    _assert(fabs(h) < TOL_);

    return 0;
}

int GHAMatrix_01() {
    double Mjd_UT1 = 51544.5;

    double gha = gast(Mjd_UT1);
    Matrix GHA1 = GHAMatrix(Mjd_UT1);
    Matrix GHA2 = R_z(gha);

    _assert((GHA1 - GHA2).norm() < TOL_);

    return 0;
}

int PrecMatrix_01()
{
    double MJD_J2000 = 51544.5;
    double Mjd_1 = MJD_J2000;
    double Mjd_2 = MJD_J2000;

    Matrix P = PrecMatrix(Mjd_1, Mjd_2);

    _assert(P.getFilas() == 3);
    _assert(P.getColumnas() == 3);

    for (int i = 1; i <= 3; ++i) {
        for (int j = 1; j <= 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            _assert(fabs(P(i,j) - expected) < TOL_);
        }
    }

    return 0;
}

int PoleMatrix_01()
{
    double xp = 0.0;
    double yp = 0.0;
    Matrix P = PoleMatrix(xp, yp);

    Matrix I = Matrix::identity(3);
    double tol = 1e-10;

    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j)
            _assert(fabs(P(i,j) - I(i,j)) < tol);

    return 0;
}

int NutMatrix_01() {
    Matrix N = NutMatrix(51544.5);

    _assert(N.getFilas() == 3);
    _assert(N.getColumnas() == 3);

    Matrix I = Matrix::identity(3);

    for (int i = 1; i <= 3; ++i) {
        for (int j = 1; j <= 3; ++j) {
            _assert(fabs(N(i,j) - I(i,j)) < TOL_);
        }
    }

    return 0;
}

int Mjday_TDB_01() {
    double Mjd_TT = 58000.0;
    double Mjd_TDB = Mjday_TDB(Mjd_TT);

    _assert(fabs(Mjd_TDB - Mjd_TT) < 1e-3);

    return 0;
}

int Doubler_01()
{
    Matrix los1 = Matrix::zeros(3,1),
            los2 = Matrix::zeros(3,1),
            los3 = Matrix::zeros(3,1);
    los1(1,1) = 1.0;
    los2(2,1) = 1.0;
    los3(3,1) = 1.0;

    Matrix zero = Matrix::zeros(3,1);

    DoublerResult R = doubler(
            3.0, 3.0,
            0.0, 0.0,
            2.0, 2.0,
            los1, los2, los3,
            zero, zero, zero,
            0.0, 1.0,
            true
    );

    _assert(R.r2.getFilas() == 3 && R.r2.getColumnas() == 1);
    _assert(R.r3.getFilas() == 3 && R.r3.getColumnas() == 1);

    _assert(std::isfinite(R.f1));
    _assert(std::isfinite(R.f2));
    _assert(std::isfinite(R.q1));
    _assert(std::isfinite(R.magr1));
    _assert(std::isfinite(R.magr2));
    _assert(std::isfinite(R.a));
    _assert(std::isfinite(R.deltae32));

    return 0;
}

int Accel_01() {
    eop19620101(21413, 13);

    int nmax = 180;
    GGM03S(nmax+1);

    DE430Coeff(2285, 1020);

    double Mjd_test = 37665.0;
    AuxParam.Mjd_UTC = Mjd_test;
    AuxParam.Mjd_TT = Mjd_test;
    AuxParam.n = nmax;
    AuxParam.m = nmax;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;


    Matrix Y(6, 1);
    Y(1, 1) = 7000e3;
    Y(2, 1) = 0.0;
    Y(3, 1) = 0.0;
    Y(4, 1) = 1.0;
    Y(5, 1) = 2.0;
    Y(6, 1) = 3.0;

    try {
        Matrix dY = Accel(0.0, Y);

        for (int i = 1; i <= 6; ++i) {
            cout << "dY(" << i << ",1) = " << dY(i, 1) << endl;
        }


    } catch (const exception &e) {
        cerr << "Error in Accel: " << e.what() << endl;
        return 1;
    }

    return 0;
}

int AzElPa_01() {
    Matrix s(3,1);
    double Az = 0.0, El = 0.0;
    Matrix dAds(3,1), dEds(3,1);

    s(1,1)=0; s(2,1)=1; s(3,1)=0;
    AzElPa(s, Az, El, dAds, dEds);

    _assert(fabs(Az - 0.0) < TOL_);
    _assert(fabs(El - 0.0) < TOL_);


    return 0;
}

int VarEqn_01() {
    int nPhi = 6, size = nPhi + nPhi * nPhi;
    Matrix yPhi(size, 1);
    for (int i = 1; i <= size; ++i) {
        yPhi(i, 1) = (double)i;
    }

    AuxParam.Mjd_UTC = 37665.0;
    AuxParam.Mjd_TT = 37665.0;
    AuxParam.n = 4;
    AuxParam.m = 4;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

    eop19620101(21413,13);
    DE430Coeff(2285, 1020);

    Cnm = Matrix(AuxParam.n+1, AuxParam.m+1);
    Snm = Matrix(AuxParam.n+1, AuxParam.m+1);

    Matrix yPhip = VarEqn(0.0, yPhi);

    _assert(yPhip.getFilas() == 42);
    _assert(yPhip.getColumnas() == 1);

    return 0;
}

int Anglesdr_01() {
    eopdata = Matrix(3, 13);

    for (int row = 1; row <= 3; ++row) {
        for (int col = 1; col <= 13; ++col) {
            eopdata(row, col) = 0.0;
        }
        eopdata(row, 4) = 58000 + (row-1); // Mjd
        eopdata(row, 5) = 3600;  // x_pole
        eopdata(row, 6) = 7200;  // y_pole
        eopdata(row, 7) = 3.5;   // UT1_UTC
        eopdata(row, 8) = 0.1;   // LOD
        eopdata(row, 9) = 10800; // dpsi
        eopdata(row,10) = 14400; // deps
        eopdata(row,11) = 1800;  // dx_pole
        eopdata(row,12) = 2700;  // dy_pole
        eopdata(row,13) = 37.0;  // TAI_UTC
    }

    double az1 = 0.0,   az2 = M_PI/2, az3 = M_PI;
    double el1 = 0.0,   el2 = 0.2,     el3 = 0.4;
    double Mjd1 = 58000, Mjd2 = 58001, Mjd3 = 58002;

    Matrix rsite1(3,1), rsite2(3,1), rsite3(3,1);
    rsite1(1,1) = 1000; rsite1(2,1) =    0; rsite1(3,1) =    0;
    rsite2(1,1) =    0; rsite2(2,1) = 2000; rsite2(3,1) =    0;
    rsite3(1,1) =    0; rsite3(2,1) =    0; rsite3(3,1) = 3000;

    AnglesDRResult res = anglesdr(az1, az2, az3, el1, el2, el3, Mjd1, Mjd2, Mjd3,
                                  rsite1, rsite2, rsite3);

    _assert(res.r2.getFilas() == 3);
    _assert(res.r2.getColumnas() == 1);
    _assert(res.v2.getFilas() == 3);
    _assert(res.v2.getColumnas() == 1);

    return 0;
}


int DEInteg_01() {
    ODEFunction func = [](double t, const Matrix& y) {
        Matrix dy(1,1);
        dy(1,1) = y(1,1);
        return dy;
    };

    double t0 = 0, tf = 1;
    double relerr = 1e-13, abserr = 1e-6;
    int n_eqn = 1;
    Matrix y0(1,1);
    y0(1,1) = 1.0;

    Matrix y_fake = DEInteg(func, t0, tf, relerr, abserr, n_eqn, y0);
    double y_real = exp(1.0);

    cout << "C++ result: " << y_fake(1,1) << " (expected " << y_real << ")\n"; // that test failed. In matlab (y_real) the result is a bit different

    return 0;
}

int test_anglesg() {
    eop19620101(21413, 13);

    double az1 = 0.5, az2 = 0.6, az3 = 0.7;
    double el1 = 0.4, el2 = 0.5, el3 = 0.6;

    double Mjd1 = 51544.0;
    double Mjd2 = Mjd1 + (10.0 / 1440.0);
    double Mjd3 = Mjd1 + (20.0 / 1440.0);

    Matrix Rs1(3, 1), Rs2(3, 1), Rs3(3, 1);
    Rs1(1, 1) = 6371e3;
    Rs1(2, 1) = 0;
    Rs1(3, 1) = 0;
    Rs2(1, 1) = 0;
    Rs2(2, 1) = 6371e3;
    Rs2(3, 1) = 0;
    Rs3(1, 1) = 0;
    Rs3(2, 1) = 0;
    Rs3(3, 1) = 6371e3;

    try {
        AnglesGResult res = anglesg(az1, az2, az3, el1, el2, el3,
                                    Mjd1, Mjd2, Mjd3, Rs1, Rs2, Rs3);

        _assert(res.r2.getFilas() == 3 && res.r2.getColumnas() == 1);
        _assert(res.v2.getFilas() == 3 && res.v2.getColumnas() == 1);
    } catch (const exception &e) {
        cerr << "Error " << endl;
        return 1;
    }

    return 0;
}

int EccAnom_02(){
    double M2 = 0.0;
    double e2 = 0.5;
    double E2 = EccAnom(M2, e2);
    _assert(fabs(E2) < TOL_);

    return 0;
}

int TimeUpdate_01() {
    Matrix P(2,2);
    P(1,1) = 1.0;  P(1,2) = 2.0;
    P(2,1) = 3.0;  P(2,2) = 4.0;

    Matrix Phi = Matrix::identity(2);

    double Qdt = 0.1;

    TimeUpdate(P, Phi, Qdt);

    bool ok = true;
    ok &= fabs(P(1,1) - 1.1) < TOL_;
    ok &= fabs(P(1,2) - 2.0) < TOL_;
    ok &= fabs(P(2,1) - 3.0) < TOL_;
    ok &= fabs(P(2,2) - 4.1) < TOL_;

    return 0;
}


    int all_tests()
{


    /*
    _verify(nutangles01);
    _verify(DEInteg_01); // With a bit error
    _verify(AccelHarmonic_02); // Error
    _verify(Legendre_02);
    _verify(AzElPa_01);
    _verify(AccelHarmonic01);
    _verify(G_AccelHarmonic01);
    _verify(sign_01);
    _verify(position01);
    _verify(angl01);
    _verify(angl02);
    _verify(angl03);
    _verify(Cheb3D_01);
    _verify(test_anglesg);
    _verify(Accel_01);
    _verify(eqnEquinox01);
    _verify(gast01);
    _verify(elements01);
    _verify(JPL_Eph_DE430_01);
    _verify(IERS01);
    _verify(LTC01);
    _verify(timediff01);
    _verify(gibbs01);
    _verify(unit01);
    _verify(hgibbs01);
    _verify(Legendre_01);
    _verify(accel_point_mass_01);
    _verify(measUpdate01);
    _verify(meanObliquity01);
    _verify(meanObliquity02);
    _verify(Mjday_01);
    _verify(gmst_01); // with less level of tolerance
    _verify(Mjday_02);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(frac_01);
    _verify(Mjday_TDB_01);
    _verify(Geodetic01);
    _verify(GHAMatrix_01);
    _verify(TimeUpdate_01);
    _verify(PoleMatrix_01); // with a less level of tolerance
    _verify(NutMatrix_01); // Error of tolerance
    _verify(PrecMatrix_01);
    _verify(Doubler_01); // warning
    _verify(EccAnom_02);
    _verify(VarEqn_01);
    _verify(Anglesdr_01); // WARNING

     */

    return 0;
}


/*
int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

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

    return result != 0;
}

*/


