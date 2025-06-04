/**
 *  @file   VarEqn.h
 *  @brief  VarEqn method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_VAREQN_H
#define PROYECTOTALLERI_VAREQN_H


#include "Matrix.h"
#include "global.h"


/**
 * @brief Computes the time derivative of the state vector and state transition matrix.
 * @param[in] x Time since epoch [s].
 * @param[in] yPhi 42×1 vector (position, velocity, 6×6 transition).
 * @return 42×1 vector.
 */
Matrix VarEqn(double x, const Matrix& yPhi);


#endif //PROYECTOTALLERI_VAREQN_H
