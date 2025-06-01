/**
 *  @file   global.h
 *  @brief  global variables
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-04-23
 ***********************************************/

#ifndef PROYECTOTALLERI_GLOBAL_H
#define PROYECTOTALLERI_GLOBAL_H

#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.h"
#include <cmath>

typedef struct{
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

extern Param AuxParam;
extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;

void eop19620101(int nfilas, int ncols);

void GGM03S(int n);

void DE430Coeff(int i, int j);

#endif


#endif //PROYECTOTALLERI_GLOBAL_H
