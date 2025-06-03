
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "Matrix.h"
#include "global.h"
#include "IERS.h"
#include "timediff.h"
#include "LTC.h"
#include "AzElPa.h"
#include "anglesg.h"
#include "MeasUpdate.h"
#include "TimeUpdate.h"
#include "gmst.h"
#include "R_z.h"
#include "Position.h"
#include "DEInteg.h"
#include "Mjday.h"
#include "Accel.h"
#include "VarEqn.h"

using namespace std;

int n_eqn;


int main() {

    try {
        eop19620101(21413, 13);
        GGM03S(180);
        DE430Coeff(2285, 1020);
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return EXIT_FAILURE;
    }


    const int nobs = 46;
    Matrix obs(nobs, 4);

    std::ifstream fobs("../data/GEOS3.txt");
    if (!fobs.is_open()) {
        std::cerr << "cant open GEOS3.txt\n";
        return EXIT_FAILURE;
    }

    auto clean  = [](const std::string& str) -> std::string {
        std::string out;
        for (char c : str)
            if (!std::isspace(static_cast<unsigned char>(c)))
                out.push_back(c);
        return out;
    };

    int i = 1;
    string line;
    while (i <= nobs && getline(fobs, line)) {

        if (line.find_first_not_of(" \t\r\n") == string::npos)
            continue;

        int  Y  = std::stoi( clean(line.substr(0,  4)) );
        int  M  = std::stoi( clean(line.substr(5,  2)) );
        int  D  = std::stoi( clean(line.substr(8,  2)) );
        int  hh = std::stoi( clean(line.substr(12, 2)) );
        int  mm = std::stoi( clean(line.substr(15, 2)) );
        double ss   = std::stod( clean(line.substr(18, 6)) );
        double az   = std::stod( clean(line.substr(25, 8)) );
        double el   = std::stod( clean(line.substr(35, 7)) );
        double Dist = std::stod( clean(line.substr(44,10)) );

        obs(i,1) = Mjday(Y, M, D, hh, mm, ss);
        obs(i,2) = Rad * az;
        obs(i,3) = Rad * el;
        obs(i,4) = 1e3 * Dist;

        ++i;
    }
    fobs.close();



    double sigma_range = 92.5;
    double sigma_az = 0.0224 * Rad;
    double sigma_el = 0.0139 * Rad;


    double lat = Rad * 21.5748;
    double lon = Rad * (-158.2706);
    double alt = 300.20;
    Matrix Rs = Position(lon, lat, alt);



    double Mjd1 = obs(1, 1), Mjd2 = obs(9, 1), Mjd3 = obs(18, 1);

    AnglesGResult gauss = anglesg(
            obs(1, 2), obs(9, 2), obs(18, 2),
            obs(1, 3), obs(9, 3), obs(18, 3),
            Mjd1, Mjd2, Mjd3, Rs, Rs, Rs
    );

    Matrix r2 = gauss.r2, v2 = gauss.v2;

    Matrix Y0_apr(6, 1);
    for (int i = 1; i <= 3; ++i) {
        Y0_apr(i, 1) = r2(i, 1);
        Y0_apr(i + 3, 1) = v2(i, 1);
    }

    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);
    double Mjd_UTC = obs(9, 1);



    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

    n_eqn = 6;



    double tLeft = -(obs(9,1) - Mjd0) * 86400.0;
    const double STEP = -600.0;

    Matrix Y = Y0_apr;
    while (std::fabs(tLeft) > 1.0) {
        double hStep = (fabs(tLeft) > fabs(STEP)) ? STEP : tLeft;
        Y = DEInteg(Accel, 0.0, hStep,
                    1e-13, 1e-6, 6, Y);
        tLeft -= hStep;
    }


    Matrix P = Matrix::zeros(6, 6);
    for (int i = 1; i <= 3; ++i) P(i, i) = 1e8;
    for (int i = 4; i <= 6; ++i) P(i, i) = 1e3;

    Matrix LT = LTC(lon, lat);


    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);


    double t = 0;
    for (int i = 1; i <= nobs; ++i) {
        double t_old = t;
        Matrix Y_old = Y;


        Mjd_UTC = obs(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;


        IERSResult ier = IERS(eopdata, Mjd_UTC, 'l');
        TimeDiffResult td = timediff(ier.UT1_UTC, ier.TAI_UTC);
        double Mjd_TT = Mjd_UTC + td.TT_UTC / 86400.0;
        double Mjd_UT1 = Mjd_TT + (ier.UT1_UTC - td.TT_UTC) / 86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;


        for (int ii = 1; ii <= 6; ++ii) {
            yPhi(ii, 1) = Y_old(ii, 1);
            for (int j = 1; j <= 6; ++j)
                yPhi(6 * j + ii, 1) = (ii == j) ? 1.0 : 0.0;
        }


        yPhi = DEInteg(VarEqn, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);


        for (int j = 1; j <= 6; ++j)
            for (int ii = 1; ii <= 6; ++ii)
                Phi(ii, j) = yPhi(6 * j + ii, 1);


        Y = DEInteg(Accel, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);


        double theta = gmst(Mjd_UT1);
        Matrix U = R_z(theta);

        Matrix r(3, 1);
        for (int i = 1; i <= 3; ++i) {
            r(i, 1) = Y(i, 1);
        }
        Matrix s = LT * (U * r - Rs);


        TimeUpdate(P, Phi);

        Matrix r3(3, 1);
        for (int ii = 1; ii <= 3; ++ii) {
            r3(ii, 1) = Y(ii, 1);
        }
        s = LT * (U * r3 - Rs);


        double Az, El;
        Matrix dAds(3, 1), dEds(3, 1);
        AzElPa(s, Az, El, dAds, dEds);


        Matrix dAds_row(1, 3);
        for (int j = 1; j <= 3; ++j) {
            dAds_row(1, j) = dAds(j, 1);
        }
        Matrix tempAz = dAds_row * LT * U;
        Matrix H_az(1, 6);
        for (int j = 1; j <= 3; ++j) {
            H_az(1, j) = tempAz(1, j);
        }
        for (int j = 4; j <= 6; ++j) {
            H_az(1, j) = 0.0;
        }


        double meas_Az = obs(i, 2);
        double pred_Az = Az;

        Matrix zAz(1,1);  zAz(1,1) = meas_Az;
        Matrix gAz(1,1);  gAz(1,1) = pred_Az;
        Matrix sAz(1,1);  sAz(1,1) = sigma_az;

        Matrix K_az = MeasUpdate(
                Y,
                zAz,
                gAz,
                sAz,
                H_az,
                P,
                6
        );


        for (int k = 1; k <= 3; ++k) r3(k,1) = Y(k,1);
        s = LT * (U * r3 - Rs);


        AzElPa(s, Az, El, dAds, dEds);


        Matrix dEdsRow(1,3);
        for (int j = 1; j <= 3; ++j) dEdsRow(1,j) = dEds(j,1);
        Matrix tempEl = dEdsRow * LT * U;

        Matrix H_el(1,6);
        for (int j = 1; j <= 3; ++j) H_el(1,j) = tempEl(1,j);
        for (int j = 4; j <= 6; ++j) H_el(1,j) = 0.0;


        Matrix zEl(1,1);  zEl(1,1) = obs(i,3);
        Matrix gEl(1,1);  gEl(1,1) = El;
        Matrix sEl(1,1);  sEl(1,1) = sigma_el;


        Matrix K_el = MeasUpdate( Y, zEl, gEl, sEl, H_el, P, 6 );





        for (int k = 1; k <= 3; ++k) r3(k,1) = Y(k,1);
        s = LT * (U * r3 - Rs);

        double Dist = sqrt( s(1,1)*s(1,1)
                            + s(2,1)*s(2,1)
                            + s(3,1)*s(3,1) );


        Matrix dDds(1,3);
        for (int j = 1; j <= 3; ++j) dDds(1,j) = s(j,1) / Dist;

        Matrix tempRg = dDds * LT * U;
        Matrix H_rg(1,6);
        for (int j = 1; j <= 3; ++j) H_rg(1,j) = tempRg(1,j);
        for (int j = 4; j <= 6; ++j) H_rg(1,j) = 0.0;


        Matrix zRg(1,1);  zRg(1,1) = obs(i,4);
        Matrix gRg(1,1);  gRg(1,1) = Dist;
        Matrix sRg(1,1);  sRg(1,1) = sigma_range;


        Matrix K_rg = MeasUpdate( Y, zRg, gRg, sRg, H_rg, P, 6 );


    }


    IERSResult ier = IERS(eopdata, obs(46, 1), 'l');
    TimeDiffResult td = timediff(ier.UT1_UTC, ier.TAI_UTC);
    double Mjd_TT = obs(46, 1) + td.TT_UTC / 86400.0;
    AuxParam.Mjd_UTC = obs(46, 1);
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix Y0 = DEInteg(Accel, 0.0, -(obs(46, 1) - obs(1, 1)) * 86400.0,
                        1e-13, 1e-6, 6, Y);

    Matrix Y_true(6, 1);
    Y_true(1, 1) = 5753.173e3;
    Y_true(2, 1) = 2673.361e3;
    Y_true(3, 1) = 3440.304e3;
    Y_true(4, 1) = 4.324207e3;
    Y_true(5, 1) = -1.924299e3;
    Y_true(6, 1) = -5.728216e3;


    cout << "\nError of Position Estimation\n";
    cout << "dX " << (Y0(1, 1) - Y_true(1, 1)) << " [m]\n";
    cout << "dY " << (Y0(2, 1) - Y_true(2, 1)) << " [m]\n";
    cout << "dZ " << (Y0(3, 1) - Y_true(3, 1)) << " [m]\n\n";
    cout << "Error of Velocity Estimation\n";
    cout << "dVx " << (Y0(4, 1) - Y_true(4, 1)) << " [m/s]\n";
    cout << "dVy " << (Y0(5, 1) - Y_true(5, 1)) << " [m/s]\n";
    cout << "dVz " << (Y0(6, 1) - Y_true(6, 1)) << " [m/s]\n";

    return 0;
}
