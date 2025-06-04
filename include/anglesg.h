/**
 *  @file   anglesg.h
 *  @brief  anglesg method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-30
 ***********************************************/

#ifndef PROYECTOTALLERI_ANGLESG_H
#define PROYECTOTALLERI_ANGLESG_H


#include "Matrix.h"

struct AnglesGResult {
    Matrix r2;
    Matrix v2;
};


/**
 * @brief this function solves the problem of orbit determination using three optical sightings.
 * @param[in] az1 azimuth at t1               rad.
 * @param[in] az2 azimuth at t2               rad.
 * @param[in] az3 azimuth at t3               rad.
 * @param[in] el1 elevation at t1             rad.
 * @param[in] el2 elevation at t2             rad.
 * @param[in] el3 elevation at t3             rad.
 * @param[in] Mjd1 Modified julian date of t1.
 * @param[in] Mjd2 Modified julian date of t2.
 * @param[in] Mjd3 Modified julian date of t3.
 * @param[in] Rs1 ijk site1 position vector   m.
 * @param[in] Rs2 ijk site2 position vector   m.
 * @param[in] Rs3 ijk site3 position vector   m.
 * @return AnglesGResult Struct containing:
 *           - r: ijk position vector at t2   m.
 *           - v: ijk velocity vector at t2   m/s.
 */
AnglesGResult anglesg(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        const Matrix& Rs1, const Matrix& Rs2, const Matrix& Rs3
);

#endif //PROYECTOTALLERI_ANGLESG_H
