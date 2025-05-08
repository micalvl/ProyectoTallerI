//
// Created by micalvl on 03/04/2025.
//

#include "TimeUpdate.h"

void TimeUpdate(Matrix& P, Matrix& Phi, double Qdt) {
    P = Phi.operator*(P)  * Phi.transpose() + Qdt;
}






