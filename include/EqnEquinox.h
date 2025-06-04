/**
 *  @file   EccAnom.h
 *  @brief  EccAnom method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#ifndef PROYECTOTALLERI_EQNEQUINOX_H
#define PROYECTOTALLERI_EQNEQUINOX_H

#include "NutAngles.h"
#include "MeanObliquity.h"
#include <cmath>


/**
 * @brief Computation of the equation of the equinoxes.
 * @param[in] Mjd_TT Modified Julian Date (Terrestrial Time).
 * @return Equation of the equinoxes.
 */
double EqnEquinox(double Mjd_TT);

#endif //PROYECTOTALLERI_EQNEQUINOX_H
