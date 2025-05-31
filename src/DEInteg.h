//
// Created by micalvl on 03/04/2025.
//

#ifndef PROYECTOTALLERI_DEINTEG_H
#define PROYECTOTALLERI_DEINTEG_H


#include "Matrix.h"
#include <functional>
#include <vector>
#include <string>

enum class DE_STATE {
    DE_INIT = 1,    // Restart integration
    DE_DONE = 2,    // Successful step
    DE_BADACC = 3,  // Accuracy requirement could not be achieved
    DE_NUMSTEPS = 4,// Permitted number of steps exceeded
    DE_STIFF = 5,   // Stiff problem suspected
    DE_INVPARAM = 6 // Invalid input parameters
};

// Alias for the function signature expected by DEInteg
using ODEFunction = std::function<Matrix(double, const Matrix&)>;


DE_STATE DEInteg(const ODEFunction& func,
                 double& t,
                 double tout,
                 double& relerr,
                 double& abserr,
                 int n_eqn,
                 Matrix& y,
                 bool PermitTOUT = true);


#endif //PROYECTOTALLERI_DEINTEG_H
