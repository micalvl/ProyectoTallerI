/**
 *  @file   Mjday_TDB.h
 *  @brief  Mjday_TDB's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-13
 ***********************************************/

#ifndef PROYECTOTALLERI_MJDAY_TDB_H
#define PROYECTOTALLERI_MJDAY_TDB_H

#include <cmath>


/**
 * @brief Computes the Modified Julian Date for barycentric dynamical time
 * @param[in] Mjd_TT Modified Julian Date in Terrestrial Time (TT).
 * @return Modified julian date (TDB)
 */
double Mjday_TDB(double Mjd_TT);


#endif //PROYECTOTALLERI_MJDAY_TDB_H
