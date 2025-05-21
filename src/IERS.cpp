//
// Created by micalvl on 03/04/2025.
//

#include "../include/IERS.h"

IERSResult IERS(const Matrix& eop, double Mjd_UTC, char interp) {
    double mjd = floor(Mjd_UTC);
    int cols = eop.getColumnas();
    int idx = 1;
    for (; idx <= cols; ++idx) {
        if (eop(4, idx) == mjd) break;
    }

    if (interp == 'l') {
        double mfme = 1440.0 * (Mjd_UTC - mjd);
        double fixf = mfme / 1440.0;
        IERSResult R;
        R.x_pole  = (eop(5,idx)  + (eop(5,idx+1)  - eop(5,idx))  * fixf) / Arcs;
        R.y_pole  = (eop(6,idx)  + (eop(6,idx+1)  - eop(6,idx))  * fixf) / Arcs;
        R.UT1_UTC =  eop(7,idx)  + (eop(7,idx+1)  - eop(7,idx))  * fixf;
        R.LOD     =  eop(8,idx)  + (eop(8,idx+1)  - eop(8,idx))  * fixf;
        R.dpsi    = (eop(9,idx)  + (eop(9,idx+1)  - eop(9,idx))  * fixf) / Arcs;
        R.deps    = (eop(10,idx) + (eop(10,idx+1) - eop(10,idx)) * fixf) / Arcs;
        R.dx_pole = (eop(11,idx) + (eop(11,idx+1) - eop(11,idx)) * fixf) / Arcs;
        R.dy_pole = (eop(12,idx) + (eop(12,idx+1) - eop(12,idx)) * fixf) / Arcs;
        R.TAI_UTC =  eop(13,idx);

        return R;
    }
    else {
        IERSResult R;
        R.x_pole  = eop(5,idx)  / Arcs;
        R.y_pole  = eop(6,idx)  / Arcs;
        R.UT1_UTC = eop(7,idx);
        R.LOD     = eop(8,idx);
        R.dpsi    = eop(9,idx)  / Arcs;
        R.deps    = eop(10,idx) / Arcs;
        R.dx_pole = eop(11,idx) / Arcs;
        R.dy_pole = eop(12,idx) / Arcs;
        R.TAI_UTC = eop(13,idx);

        return R;
    }
}
