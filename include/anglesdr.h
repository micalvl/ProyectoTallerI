/**
 *  @file   anglesdr.h
 *  @brief  anglesdr function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-30
 ***********************************************/

#ifndef PROYECTOTALLERI_ANGLESDR_H
#define PROYECTOTALLERI_ANGLESDR_H

#include "Matrix.h"

struct AnglesDRResult {
    Matrix r2;
    Matrix v2;
};


/**
 * @brief Determines geocentric position and velocity from three line-of-sight observations using the double-r iteration method.
 *
 * Given three azimuth/elevation measurements at distinct epochs and the corresponding site position vectors,
 * this function iteratively refines the ranges (double-r) to compute the spacecraft’s geocentric position and velocity
 * at the time of the second observation.
 *
 * @param[in] az1    Azimuth angle of first observation [rad].
 * @param[in] az2    Azimuth angle of second observation [rad].
 * @param[in] az3    Azimuth angle of third observation [rad].
 * @param[in] el1    Elevation angle of first observation [rad].
 * @param[in] el2    Elevation angle of second observation [rad].
 * @param[in] el3    Elevation angle of third observation [rad].
 * @param[in] Mjd1   Modified Julian Date of first observation (UTC).
 * @param[in] Mjd2   Modified Julian Date of second observation (UTC).
 * @param[in] Mjd3   Modified Julian Date of third observation (UTC).
 * @param[in] rsite1 3×1 geocentric site position vector at time Mjd1.
 * @param[in] rsite2 3×1 geocentric site position vector at time Mjd2.
 * @param[in] rsite3 3×1 geocentric site position vector at time Mjd3.
 * @return AnglesDRResult
 *         Structure containing:
 *           - r2: 3×1 spacecraft position vector at time Mjd2 (geocentric coordinates).
 *           - v2: 3×1 spacecraft velocity vector at time Mjd2 (geocentric coordinates).
 */
AnglesDRResult anglesdr(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        const Matrix& rsite1,
        const Matrix& rsite2,
        const Matrix& rsite3
);


#endif //PROYECTOTALLERI_ANGLESDR_H
