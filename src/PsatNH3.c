// PsatNH3.c

#include "lagrangian/env.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <PsatNH3.h>

/// @brief Calculates Saturaiton Pressure for Ammonia based on Temperature, based on Data from NIST webbook (Stull 1947)
/// @param Td - Temperatuer of droplet
/// @param P_sat - Saturation Pressure
void Saturation_PressureNH3(CONVERGE_precision_t Td, CONVERGE_precision_t *P_sat)
{
  CONVERGE_precision_t psatA1, psatB1, psatC1, psatA2, psatB2, psatC2;
  CONVERGE_precision_t Lim1, Lim2, Lim3;
  Lim3 = 371.5;
  Lim2 = 239.6;
  Lim1 = 164.0;
  //Higher temp range constants
  psatA1 = 3.18757;
  psatB1 = 506.713;
  psatC1 = -80.78;
  //lower temp range constants 
  psatA2 = 4.86886;
  psatB2 = 1113.928;
  psatC2 = -10.409;
  if ((Td > Lim2) & (Td < Lim3))
  {
    *P_sat = pow(10, (psatA2 - (psatB2 / (psatC2 + Td)))) * 1e5; // P_sat in Pa	from NIST Webbook for iso-octane
  }
  else if ((Td < Lim2) & (Td > Lim1))
  {
     *P_sat = pow(10, (psatA1 - (psatB1 / (psatC1 + Td)))) * 1e5; // P_sat in Pa	from NIST Webbook for iso-octane
  }
  else if (Td > Lim3)
  {
    printf("\nTd too hot for correlation (T=%f K), aborting....", Td);
    CONVERGE_mpi_abort();
  }
  else if (Td < Lim1)
  {
    printf("\n Td too cold for correlation (T=%f K), aborting...", Td);
    CONVERGE_mpi_abort();
  }
}