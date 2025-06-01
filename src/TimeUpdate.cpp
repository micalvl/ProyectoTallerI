/**
 *  @file   TimeUpdate.cpp
 *  @brief  TimeUpdate's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-14
 ***********************************************/

#include "../include/TimeUpdate.h"

void TimeUpdate(Matrix& P, const Matrix& Phi, double Qdt) {
    Matrix Q = Matrix::identity(P.getFilas()).opsc(Qdt);
    P = Phi.operator*(P)  * Phi.transpose() + Q;
}






