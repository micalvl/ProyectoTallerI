/**
 *  @file   TimeUpdate.h
 *  @brief  TimeUpdate's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-14
 ***********************************************/

#ifndef PROYECTOTALLERI_TIMEUPDATE_H
#define PROYECTOTALLERI_TIMEUPDATE_H

#include "Matrix.h"

/**
 * @brief Performs the time update (prediction) of a covariance matrix.
 * @param[in,out] P Current covariance matrix (n×n). On return, P is replaced by the updated covariance.
 * @param[in] Phi State transition matrix (n×n).
 * @param[in] Qdt Process noise covariance increment (scalar), added as Q = I·Qdt.
 */
void TimeUpdate(Matrix& P, const Matrix& Phi, double Qdt = 0.0);


#endif //PROYECTOTALLERI_TIMEUPDATE_H
