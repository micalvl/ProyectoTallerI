/**
 *  @file   unit.cpp
 *  @brief  unit method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-14
 ***********************************************/

#include "../include/unit.h"
#include <cmath>

Matrix unit(const Matrix& vec) {
    const double small = 0.000001;
    double magv = vec.norm();
    Matrix out = Matrix::zeros(vec.getFilas(), vec.getColumnas());
    if (magv > small) {
        for (int i = 1; i <= vec.getFilas(); ++i) {
            for (int j = 1; j <= vec.getColumnas(); ++j) {
                out(i,j) = vec(i,j) / magv;
            }
        }
    }
    return out;
}