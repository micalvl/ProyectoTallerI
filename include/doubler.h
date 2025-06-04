/**
 *  @file   doubler.h
 *  @brief  doubler method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_DOUBLER_H
#define PROYECTOTALLERI_DOUBLER_H

#include "Matrix.h"
#include "Sat_const.h"


struct DoublerResult {
    Matrix r2, r3;
    double f1, f2, q1;
    double magr1, magr2, a, deltae32;
};


/**
 * @brief Computes the position vectors and orbital parameters using three site observations.
 * @param[in] cc1 Coefficient associated with first observation.
 * @param[in] cc2 Coefficient associated with second observation.
 * @param[in] magrsite1 Magnitude of first site position vector.
 * @param[in] magrsite2 Magnitude of second site position vector.
 * @param[in] magr1in Initial magnitude guess for first position vector.
 * @param[in] magr2in Initial magnitude guess for second position vector.
 * @param[in] los1 Line-of-sight unit vector for the first observation.
 * @param[in] los2 Line-of-sight unit vector for the second observation.
 * @param[in] los3 Line-of-sight unit vector for the third observation.
 * @param[in] rsite1 Site position vector for first observation.
 * @param[in] rsite2 Site position vector for second observation.
 * @param[in] rsite3 Site position vector for third observation.
 * @param[in] t1 Time associated with first observation.
 * @param[in] t3 Time associated with third observation.
 * @param[in] direct Boolean flag indicating the direction of motion (true if direct).
 * @return A DoublerResult struct
 */
DoublerResult doubler(
        double cc1, double cc2,
        double magrsite1, double magrsite2,
        double magr1in,   double magr2in,
        const Matrix& los1, const Matrix& los2, const Matrix& los3,
        const Matrix& rsite1,const Matrix& rsite2,const Matrix& rsite3,
        double t1, double t3,
        bool direct
);


#endif //PROYECTOTALLERI_DOUBLER_H
