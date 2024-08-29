//DGRE_NH3.h

#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>


#ifndef DGRE_NH3_H
#define DGRE_NH3_H


void DGRE_NH3(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx,CONVERGE_precision_t global_density);


#endif // DGRE_NH3_H
