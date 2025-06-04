/**
 *  @file   timediff.cpp
 *  @brief  timediff method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/


#include "../include/TimeDiff.h"

TimeDiffResult timediff(double UT1_UTC, double TAI_UTC) {
    const double TT_TAI  = +32.184;   // TT−TAI difference[s]
    const double GPS_TAI = -19.0;     // GPS−TAI difference [s]

    double TT_GPS  = TT_TAI - GPS_TAI;  // TT−GPS difference
    double TAI_GPS = -GPS_TAI;          // TAI−GPS difference
    double UT1_TAI = UT1_UTC - TAI_UTC; // UT1−TAI difference
    double UTC_TAI = -TAI_UTC;          // UTC−TAI difference
    double UTC_GPS = UTC_TAI - GPS_TAI; // UTC−GPS difference
    double UT1_GPS = UT1_TAI - GPS_TAI; // UT1−GPS difference
    double TT_UTC  = TT_TAI - UTC_TAI;  // TT−UTC difference
    double GPS_UTC = GPS_TAI - UTC_TAI; // GPS−UTC difference

    return { UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC };
}