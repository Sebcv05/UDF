//RPE_euler.h
// Rayleigh-Plesset Equation solver with Explicit Euler integration
// For bubble growth within finite liquid droplets

#ifndef RPE_EULER_H
#define RPE_EULER_H

#include "lagrangian/env.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// Physical properties structure for RPE solver
typedef struct {
    // Liquid properties
    CONVERGE_precision_t rho_l;       // Liquid density (kg/m³)
    CONVERGE_precision_t mu_l;        // Dynamic viscosity (Pa·s)
    CONVERGE_precision_t k_l;         // Thermal conductivity (W/(m·K))
    CONVERGE_precision_t cp_l;        // Specific heat (J/(kg·K))
    CONVERGE_precision_t sigma;       // Surface tension (N/m)
    
    // Vapor properties
    CONVERGE_precision_t R_spec;      // Specific gas constant (J/(kg·K))
    CONVERGE_precision_t L_v;         // Latent heat (J/kg)
    
    // Geometry
    CONVERGE_precision_t Ro;          // Droplet radius (m)
    CONVERGE_precision_t m_drop;      // Droplet mass (kg)
    
    // Environment
    CONVERGE_precision_t P_amb;       // Ambient pressure (Pa)
    CONVERGE_precision_t T_amb;       // Ambient temperature (K)
    
    // Numerical limits
    CONVERGE_precision_t max_Nu;      // Maximum Nusselt number
} RPE_Params;

// Bubble state structure
typedef struct {
    CONVERGE_precision_t R;           // Bubble radius (m)
    CONVERGE_precision_t Rdot;        // Bubble wall velocity (m/s)
    CONVERGE_precision_t T_drop;      // Droplet temperature (K)
    CONVERGE_precision_t m_b;         // Bubble vapor mass (kg)
} BubbleState;

// Derivative structure
typedef struct {
    CONVERGE_precision_t dRdt;
    CONVERGE_precision_t dRdotdt;
    CONVERGE_precision_t dTdt;
    CONVERGE_precision_t dmbdt;
} BubbleDerivatives;

// Main solver function - replaces Bubble_Velocity()
void RPE_euler_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
);

// Core computation functions
void compute_thermal_mass_transfer(
    CONVERGE_precision_t R,
    CONVERGE_precision_t Rdot,
    CONVERGE_precision_t T_drop,
    CONVERGE_precision_t m_b,
    const RPE_Params* params,
    CONVERGE_precision_t* Nu_out,
    CONVERGE_precision_t* Q_out,
    CONVERGE_precision_t* mdot_out
);

void compute_derivatives(
    const BubbleState* state,
    const RPE_Params* params,
    BubbleDerivatives* derivs,
    CONVERGE_precision_t mdot
);

void euler_step(
    BubbleState* state,
    const BubbleDerivatives* derivs,
    CONVERGE_precision_t dt
);

// Utility functions
CONVERGE_precision_t safe_sqrt(CONVERGE_precision_t x);
CONVERGE_precision_t safe_divide(CONVERGE_precision_t num, CONVERGE_precision_t denom, CONVERGE_precision_t fallback);

#endif // RPE_EULER_H
