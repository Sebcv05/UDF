// Breakup_Song.h
// Header for Song et al. breakup model

#ifndef BREAKUP_SONG_H
#define BREAKUP_SONG_H

#include "lagrangian/env.h"
#include <CONVERGE/udf.h>

// ============================================================================
// Song Breakup Function Declaration
// ============================================================================
// Performs breakup for Song et al. isothermal bubble growth model
// - Randomly creates 2-5 child droplets (equal probability)
// - R_child = R_parent / cbrt(N_child_droplets)
// - num_drop_child = num_drop_parent × N_child_droplets (volume balance)
// - Updates parent parcel in-place (no new parcel creation)
// - Calculates radial velocity from momentum balance
// ============================================================================
void Breakup_Song(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb
);

#endif // BREAKUP_SONG_H
