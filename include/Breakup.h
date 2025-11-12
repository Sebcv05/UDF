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

// Rosin-Rammler distribution parameters
typedef struct {
    double n_RR;              // Shape parameter (2.5-4.0)
    double gamma_ratio;       // Pre-computed tgamma(1+2/n) / tgamma(1+3/n)
    int initialized;          // Flag to ensure one-time initialization
} RR_Params;

void init_RR_distribution(double n_shape);
void Breakup(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx,CONVERGE_cloud_t cloud);

#endif // BREAKUP_H