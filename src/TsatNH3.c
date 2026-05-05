//TsatNH3.c

#include "lagrangian/env.h"
#include <user_header.h>
#include <spray_break.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>


/// @brief Calculates Saturaiton Temperatue of Ammonia given ambient Pressure (NIST Webbook, Stull 1947)
/// @param P - Ambient Pressure
/// @return - Saturation Temperature 
CONVERGE_precision_t T_satNH3(CONVERGE_precision_t P)
   {



      CONVERGE_precision_t tsatA1,tsatA2,tsatB1, tsatB2,tsatC1, tsatC2, P_bar;
      tsatA1 = 4.86886;
      tsatB1 = 1113.928;
      tsatC1 = -10.409;
      tsatA2 = 3.18757;
      tsatB2 = 506.713;
      tsatC2 = -80.78;
      P_bar = P / 1e5;
      if(P_bar> 0.9933080)
      {
        return (tsatB1 / (tsatA1 - log10(P_bar))) - tsatC1;
      }else
      {
        return (tsatB2 / (tsatA2 - log10(P_bar))) - tsatC2;
      }
      
      // printf("Tsat = %e\n", (B / (A - log10(P_bar))) - C);
   }