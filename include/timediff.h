/**
 *  @file   timediff.h
 *  @brief  timediff method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_TIMEDIFF_H
#define PROYECTOTALLERI_TIMEDIFF_H


struct TimeDiffResult {
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
};

/**
 * @brief Computes standard time differences between time scales.
 * @param[in] UT1_UTC Difference UT1−UTC [s].
 * @param[in] TAI_UTC Difference TAI−UTC [s].
 * @return A TimeDiffResult struct
 */
TimeDiffResult timediff(double UT1_UTC, double TAI_UTC);


#endif //PROYECTOTALLERI_TIMEDIFF_H
