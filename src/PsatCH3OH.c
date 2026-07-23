// PsatNH3.c

#include "lagrangian/env.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <PsatNH3.h>

/// @brief Calculates Saturaiton Pressure for Methanol based on Temperature, based on Data from NIST webbook (Stull 1947)
/// @param Td - Temperatuer of droplet
/// @param P_sat - Saturation Pressure
void Saturation_PressureCH3OH(CONVERGE_precision_t Td, CONVERGE_precision_t *P_sat)
{
  CONVERGE_precision_t psatA1, psatB1, psatC1, psatA2, psatB2, psatC2;
  CONVERGE_precision_t Lim1, Lim2, Lim3;
  Lim3 = 512.63;
  Lim2 = 353.5;
  Lim1 = 288.1;
  //Higher temp range constants
  psatA1 = 5.15853;
  psatB1 = 1569.613;
  psatC1 = -34.846;
  //lower temp range constants 
  psatA2 = 5.20509;
  psatB2 = 1581.341;
  psatC2 = -33.5;
  if ((Td > Lim2) && (Td < Lim3))
  {
    *P_sat = pow(10, (psatA2 - (psatB2 / (psatC2 + Td)))) * 1e5; // P_sat in Pa	from NIST Webbook for iso-octane
  }
  else if ((Td < Lim2) && (Td > Lim1))
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