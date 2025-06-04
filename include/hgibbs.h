/**
 *  @file   hgibbs.h
 *  @brief  hgibbs method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-17
 ***********************************************/

#ifndef PROYECTOTALLERI_HGIBBS_H
#define PROYECTOTALLERI_HGIBBS_H


#include "Matrix.h"
#include <string>
#include "Sat_const.h"
#include "angl.h"
#include "unit.h"

struct HGibbsResult {
    Matrix    v2;
    double    theta;
    double    theta1;
    double    copa;
    std::string error;
};


/**
 * @brief This function implements the Herrick–Gibbs approximation for orbit determination, and finds the middle velocity vector for the three given position vectors.
 * @param[in] r1 ijk position vector #1 [m]
 * @param[in] r2 ijk position vector #2 [m]
 * @param[in] r3 ijk position vector #3 [m]
 * @param[in] Mjd1 Julian date of 1st sighting [days from 4713 BC]
 * @param[in] Mjd2 Julian date of 2nd sighting [days from 4713 BC]
 * @param[in] Mjd3 Julian date of 3rd sighting [days from 4713 BC]
 * @return  HGibbsResult struct
 */
HGibbsResult hgibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3, double Mjd1, double Mjd2, double Mjd3);

#endif //PROYECTOTALLERI_HGIBBS_H
