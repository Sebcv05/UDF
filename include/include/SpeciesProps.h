//SpeciesProps.h

#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>


#ifndef SPECIES_PROPS_H
#define SPECIES_PROPS_H


void species_props(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx, CONVERGE_size_t num_parcel_species, CONVERGE_precision_t* H, CONVERGE_precision_t* csubp_l,CONVERGE_species_t sp);


#endif // SPECIES_PROPS_H