/**
 *  @file   sign_.cpp
 *  @brief  sign_'s method
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-05-08
 ***********************************************/

#include <cmath>
#include "../include/sign_.h"


double sign_(double a, double b){
    if (b>=0.0){
        return fabs(a);
    }
    else{
        return - fabs(a);
    }

}


