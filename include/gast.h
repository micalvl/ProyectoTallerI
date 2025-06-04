/**
 *  @file   gast.h
 *  @brief  gast method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/
#ifndef PROYECTOTALLERI_GAST_H
#define PROYECTOTALLERI_GAST_H

#include "gmst.h"
#include "EqnEquinox.h"
#include <cmath>
using namespace std;


/**
 * @brief Greenwich Apparent Sidereal Time.
 * @param[in] Mjd_UT1 Modified Julian Date UT1.
 * @return GAST in [rad]
 */
double gast(double Mjd_UT1);


#endif //PROYECTOTALLERI_GAST_H
