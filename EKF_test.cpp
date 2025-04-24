#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include "../include/Mjday.h"
#include "../include/Matrix.h"
#include "../include/R_x.h"
#include "../include/Legendre.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/gmst.h"
#include "../include/Frac.h"

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

int Legendre_01(){
    Matrix pnm(3,3), dpnm(3,3);
    Legendre(2,2,1.0, pnm, dpnm);
    pnm.print();
    dpnm.print();
    return 0;
}

double gmst_01() {

    double Mjd_24_04_2025 = 60825.0;

    double gmst_rad = gmst(Mjd_24_04_2025);

    double expected_gmst = 3.8397;

    _assert(fabs(gmst_rad - expected_gmst) < TOL_);

    return 0.0;
}

int frac_01(){
    _assert(fabs(Frac(1.5) - 0.5) < TOL_);
    _assert(fabs(Frac(2.25) - 0.25) < TOL_);
    _assert(fabs(Frac(1.66) - 0.66) < TOL_);
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



