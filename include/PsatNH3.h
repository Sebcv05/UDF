// PsatNH3.h

#include "lagrangian/env.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>


#ifndef P_SAT_NH3_H
#define P_SAT_NH3_H

void Saturation_PressureNH3(CONVERGE_precision_t Td,CONVERGE_precision_t* P_sat);

#endif // P_SAT_NH3_H