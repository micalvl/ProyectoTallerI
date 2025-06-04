/**
 *  @file   EccAnom.h
 *  @brief  EccAnom method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-02
 ***********************************************/

#ifndef PROYECTOTALLERI_ECCANOM_H
#define PROYECTOTALLERI_ECCANOM_H


/**
 * @brief Computes the eccentric anomaly for elliptic orbits.
 * @param[in] M Mean anomaly in [rad].
 * @param[in] e Eccentricity of the orbit [0,1].
 * @return Eccentric anomaly in [rad].
 */
double EccAnom(double M, double e);

#endif //PROYECTOTALLERI_ECCANOM_H
