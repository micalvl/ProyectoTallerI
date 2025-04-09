/**
 *  @file   Mjday.h
 *  @brief  Mjday's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-03
 ***********************************************/

#ifndef PROYECTOTALLERI_MJDAY_H
#define PROYECTOTALLERI_MJDAY_H

/**
 * @brief Calculates the Modified Julian Date (MJD) from a calendar date and time.
 *
 * @param yr The year of the date.
 * @param mon The month of the date.
 * @param day The day of the date.
 * @param hr The hour of the time (default 0).
 * @param min The minute of the time (default 0).
 * @param sec The second of the time (default 0).
 * @return The Modified Julian Date.
 */
double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0);


#endif //PROYECTOTALLERI_MJDAY_H
