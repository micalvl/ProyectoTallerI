/**
 *  @file   G_AccelHarmonic.h
 *  @brief  G_AccelHarmonic function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-07
 ***********************************************/

#ifndef PROYECTOTALLERI_G_ACCELHARMONIC_H
#define PROYECTOTALLERI_G_ACCELHARMONIC_H


#include "Matrix.h"
#include "AccelHarmonic.h"
/**
 * @brief Calcula el gradiente del campo gravitatorio armónico de la Tierra.
 *
 * @param r Vector de posición del satélite en el sistema verdadero de fecha.
 * @param U Matriz de transformación al sistema cuerpo-fijo.
 * @param n_max Grado máximo del modelo gravitatorio.
 * @param m_max Orden máximo del modelo gravitatorio.
 * @return Matriz 3x3 con el gradiente (G = da/dr) en el sistema verdadero de fecha.
 */
Matrix G_AccelHarmonic(const Matrix& r, const Matrix& E, int n_max, int m_max);


#endif //PROYECTOTALLERI_G_ACCELHARMONIC_H
