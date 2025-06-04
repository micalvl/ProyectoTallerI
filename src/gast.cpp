/**
 *  @file   gast.cpp
 *  @brief  gast method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   undefined
 ***********************************************/

#include "gast.h"


double gast (double Mjd_UT1){
    double gstime = fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*M_PI );
    return gstime;
}



