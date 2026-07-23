//BubbleDensityCH3OH.c

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
#include<BubbleDensityCH3OH.h>

   CONVERGE_precision_t bubble_densityCH3OH(CONVERGE_precision_t P, CONVERGE_precision_t T)
   {

      CONVERGE_precision_t R = 8.3145 * 1000 / 32.0419 ; // R for CH3OH from NIST webbook
      CONVERGE_precision_t rho_b = P / (R * T);          // Assuming Ideal Gas
      return (rho_b);
   }