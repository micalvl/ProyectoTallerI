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

void TimeUpdate(Matrix& P, const Matrix& Phi, double Qdt = 0.0);


#endif //PROYECTOTALLERI_TIMEUPDATE_H
