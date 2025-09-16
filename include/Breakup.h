//Breakup.h

#include "lagrangian/env.h"

#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <globals.h>
#include "CONVERGE/tools/vec3.h"

#ifndef BREAKUP_H
#define BREAKUP_H

void Breakup(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx,CONVERGE_cloud_t cloud);

#endif // BREAKUP_H