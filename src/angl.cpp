/**
 *  @file   angl.cpp
 *  @brief  angl method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-13
 ***********************************************/

#include "../include/angl.h"

double angl(const Matrix& vec1, const Matrix& vec2) {
    const double small     = 1e-8;
    const double undefined = 999999.1;

    double magv1 = vec1.norm();
    double magv2 = vec2.norm();

    if (magv1 * magv2 > small * small) {
        double dot = 0.0;
        int n = vec1.getFilas();
        for (int i = 1; i <= n; ++i) {
            dot += vec1(i,1) * vec2(i,1);
        }
        double temp = dot / (magv1 * magv2);
        if (fabs(temp) > 1.0) {
            temp = (temp > 0.0 ? 1.0 : -1.0);
        }
        return acos(temp);
    } else {
        return undefined;
    }
}