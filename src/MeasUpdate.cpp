//
// Created by micalvl on 03/04/2025.
//

#include "MeasUpdate.h"


Matrix MeasUpdate(
        Matrix&       x,
        const Matrix& z,
        const Matrix& g,
        const Matrix& s,
        const Matrix& G,
        Matrix&       P,
        int           n
) {
    int m = z.getFilas();

    Matrix Inv_W(m, m);
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= m; ++j)
            Inv_W(i,j) = 0.0;
        Inv_W(i,i) = s(i,1) * s(i,1);
    }

    Matrix GT   = G.transpose();
    Matrix S    = Inv_W + G * P * GT;
    Matrix Sinv = S.inverse();
    Matrix K    = P * GT * Sinv;

    x = x + K * (z - g);
    Matrix I    = Matrix::identity(n);
    P = (I - K * G) * P;

    return K;
}