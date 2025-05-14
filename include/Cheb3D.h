/**
 *  @file   Cheb3D.h
 *  @brief  Cheb's function
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-06
 ***********************************************/

#ifndef PROYECTOTALLERI_CHEB3D_H
#define PROYECTOTALLERI_CHEB3D_H

#include <iostream>
#include <vector>
using namespace std;

vector<double> Cheb3D(double t, int N, double Ta, double Tb,const vector<double>& Cx,const vector<double>& Cy,const vector<double>& Cz);


#endif //PROYECTOTALLERI_CHEB3D_H
