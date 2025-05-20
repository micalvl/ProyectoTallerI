//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_TIMEDIFF_H
#define PROYECTOTALLERI_TIMEDIFF_H


struct TimeDiffResult {
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
};

TimeDiffResult timediff(double UT1_UTC, double TAI_UTC);


#endif //PROYECTOTALLERI_TIMEDIFF_H
