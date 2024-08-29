//Vb.h

#include "lagrangian/env.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>


#ifndef V_B_H
#define V_B_H

void Bubble_Velocity(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx,CONVERGE_precision_t P_sat, CONVERGE_precision_t P_amb);

#endif // V_B_H