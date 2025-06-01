//
// Created by micalvl on 03/04/2025.
//

#include <stdexcept>
#include <iostream>
#include "../include/IERS.h"
using namespace std;

IERSResult IERS(const Matrix& eop, double Mjd_UTC, char interp) {
    // Verificación mínima necesaria de dimensiones
    if (eop.getColumnas() < 13) {
        throw std::runtime_error("IERS: eopdata have min 13 columns");
    }

    double mjd = std::floor(Mjd_UTC);
    int rows = eop.getFilas();
    int idx = 1;

    for (; idx <= rows; ++idx) {
        if (fabs(eop(idx, 4) - mjd) < 1e-8) break;
    }

    if (idx > rows) {
        throw std::runtime_error("IERS: no existe fila con floor(Mjd_UTC) en eopdata.");
    }

    IERSResult R;
    std::cout << "[DEBUG] idx = " << idx << ", rows = " << rows << ", interp = " << interp << std::endl;

    if (interp == 'l' && idx < rows) {
        std::cout << "[DEBUG] Acceso interpolación: idx+1 = " << (idx+1) << std::endl;
        double mfme = 1440.0 * (Mjd_UTC - mjd);
        double fixf = mfme / 1440.0;

        R.x_pole  = (eop(idx, 5)   + (eop(idx+1, 5)  - eop(idx, 5))  * fixf) / Arcs;
        R.y_pole  = (eop(idx, 6)   + (eop(idx+1, 6)  - eop(idx, 6))  * fixf) / Arcs;
        R.UT1_UTC = eop(idx, 7)   + (eop(idx+1, 7)  - eop(idx, 7))  * fixf;
        R.LOD     = eop(idx, 8)   + (eop(idx+1, 8)  - eop(idx, 8))  * fixf;
        R.dpsi    = (eop(idx, 9)   + (eop(idx+1, 9)  - eop(idx, 9))  * fixf) / Arcs;
        R.deps    = (eop(idx,10)   + (eop(idx+1,10)  - eop(idx,10))  * fixf) / Arcs;
        R.dx_pole = (eop(idx,11)   + (eop(idx+1,11)  - eop(idx,11))  * fixf) / Arcs;
        R.dy_pole = (eop(idx,12)   + (eop(idx+1,12)  - eop(idx,12))  * fixf) / Arcs;
        R.TAI_UTC = eop(idx,13);
        return R;
    }

    R.x_pole  = eop(idx, 5) / Arcs;
    R.y_pole  = eop(idx, 6) / Arcs;
    R.UT1_UTC = eop(idx, 7);
    R.LOD     = eop(idx, 8);
    R.dpsi    = eop(idx, 9) / Arcs;
    R.deps    = eop(idx,10) / Arcs;
    R.dx_pole = eop(idx,11) / Arcs;
    R.dy_pole = eop(idx,12) / Arcs;
    R.TAI_UTC = eop(idx,13);
    return R;
}
