//RPE_song.h
// Song et al. Rayleigh-Plesset Equation solver (isothermal model)

#ifndef RPE_SONG_H
#define RPE_SONG_H

#include "lagrangian/env.h"
#include <CONVERGE/udf.h>

// Song RPE model parameters
typedef struct {
    CONVERGE_precision_t rho_l;       // Liquid density (kg/m³)
    CONVERGE_precision_t mu_l;        // Dynamic viscosity (Pa·s)
    CONVERGE_precision_t sigma;       // Surface tension (N/m)
    CONVERGE_precision_t R_spec;      // Specific gas constant (488.2 J/(kg·K) for NH3)
    CONVERGE_precision_t P_amb;       // Ambient pressure (Pa)
    CONVERGE_precision_t P_r0;        // Residual gas pressure (Pa) - typically 1e6
    CONVERGE_precision_t kappa;       // Surface viscosity (Pa·s·m) - typically 0.0
} SongParams;

// Main solver function
void RPE_song_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
);

// Helper functions
CONVERGE_precision_t song_compute_void_fraction(
    CONVERGE_precision_t R_bubble,
    CONVERGE_precision_t R_drop_current
);

CONVERGE_precision_t song_compute_mixture_density(
    CONVERGE_precision_t epsilon,
    CONVERGE_precision_t rho_v,
    CONVERGE_precision_t rho_l
);

CONVERGE_precision_t song_compute_acceleration(
    CONVERGE_precision_t R,
    CONVERGE_precision_t Rdot,
    CONVERGE_precision_t R0,
    CONVERGE_precision_t rho_m,
    CONVERGE_precision_t P_sat,
    const SongParams* params
);

#endif // RPE_SONG_H
