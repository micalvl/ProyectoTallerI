/**
 *  @file   Accel.h
 *  @brief  Accel's method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-30
 ***********************************************/

#ifndef PROYECTOTALLERI_ACCEL_H
#define PROYECTOTALLERI_ACCEL_H


#include "global.h"
#include "Matrix.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Mjday_TDB.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "Sat_const.h"

using namespace std;

Matrix Accel(double x, const Matrix& Y);

#endif //PROYECTOTALLERI_ACCEL_H
