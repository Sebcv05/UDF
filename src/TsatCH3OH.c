//TsatCH3OH.c

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

/// @brief Calculates Saturaiton Temperatue of Methanol given ambient Pressure (NIST Webbook, Stull 1947)
/// @param P - Ambient Pressure
/// @return - Saturation Temperature 
CONVERGE_precision_t T_satCH3OH(CONVERGE_precision_t P)
   {



      CONVERGE_precision_t tsatA1,tsatA2,tsatB1, tsatB2,tsatC1, tsatC2, P_bar;
      tsatA1 = 5.15853;
      tsatB1 = 1569.613;
      tsatC1 = -34.846;
      tsatA2 = 5.20509;
      tsatB2 = 1581.341;
      tsatC2 = -33.5;
      P_bar = P / 1e5;
      if(P_bar> 1.70911)
      {
        return (tsatB1 / (tsatA1 - log10(P_bar))) - tsatC1;
      }else
      {
        return (tsatB2 / (tsatA2 - log10(P_bar))) - tsatC2;
      }
      
      // printf("Tsat = %e\n", (B / (A - log10(P_bar))) - C);
   }