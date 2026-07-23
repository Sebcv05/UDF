//DGRE_CH3OH.h

#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>


#ifndef DGRE_CH3OH_H
#define DGRE_CH3OH_H


void DGRE_CH3OH(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx,CONVERGE_precision_t global_density);


#endif // DGRE_CH3OH_H
