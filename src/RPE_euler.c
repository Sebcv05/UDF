//RPE_euler.c
// Rayleigh-Plesset Equation solver with Explicit Euler integration

#include "lagrangian/env.h"
#include <RPE_euler.h>
#include <PsatNH3.h>
#include <TsatNH3.h>
#include <BubbleDensityNH3.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>

// Physical constants
#define PI 3.14159265358979323846
#define R_SPEC_NH3 488.2  // J/(kg·K) - Specific gas constant for ammonia

// Utility functions
CONVERGE_precision_t safe_sqrt(CONVERGE_precision_t x) {
    return (x > 0.0) ? CONVERGE_sqrt(x) : 0.0;
}

CONVERGE_precision_t safe_divide(CONVERGE_precision_t num, CONVERGE_precision_t denom, CONVERGE_precision_t fallback) {
    return (fabs(denom) > 1e-20) ? (num / denom) : fallback;
}

// Compute heat and mass transfer rates
void compute_thermal_mass_transfer(
    CONVERGE_precision_t R,
    CONVERGE_precision_t Rdot,
    CONVERGE_precision_t T_drop,
    CONVERGE_precision_t m_b,
    const RPE_Params* params,
    CONVERGE_precision_t* Nu_out,
    CONVERGE_precision_t* Q_out,
    CONVERGE_precision_t* mdot_out
) {
    // Calculate bubble pressure from current mass
    CONVERGE_precision_t Vb = (4.0/3.0) * PI * R*R*R;
    CONVERGE_precision_t rho_v = safe_divide(m_b, Vb, 1e-6);
    if (rho_v < 1e-6) rho_v = 1e-6;
    
    CONVERGE_precision_t T_drop_safe = (T_drop > 1e-3) ? T_drop : 1e-3;
    CONVERGE_precision_t Pb = rho_v * params->R_spec * T_drop_safe;
    
    // Interface temperature from bubble pressure
    CONVERGE_precision_t T_int = T_satNH3(Pb);
    
    // Temperature difference for heat transfer
    CONVERGE_precision_t dT = T_drop - T_int;
    
    // Reynolds and Prandtl numbers
    CONVERGE_precision_t L_char = 2.0 * R;
    CONVERGE_precision_t Re = params->rho_l * fabs(Rdot) * L_char / params->mu_l;
    CONVERGE_precision_t Pr = params->mu_l * params->cp_l / params->k_l;
    
    // Nusselt number (Ranz-Marshall correlation)
    CONVERGE_precision_t Nu_uncapped = 2.0 + 0.6 * safe_sqrt(Re) * pow(Pr, 1.0/3.0);
    CONVERGE_precision_t Nu = (Nu_uncapped < params->max_Nu) ? Nu_uncapped : params->max_Nu;
    
    // Heat transfer coefficient and rate
    CONVERGE_precision_t h_conv = Nu * params->k_l / L_char;
    CONVERGE_precision_t A_bubble = 4.0 * PI * R * R;
    CONVERGE_precision_t Q_conv = h_conv * A_bubble * dT;
    
    // Mass transfer rate (thermal limiting)
    CONVERGE_precision_t mdot = safe_divide(Q_conv, params->L_v, 0.0);
    
    // Output
    *Nu_out = Nu;
    *Q_out = Q_conv;
    *mdot_out = mdot;
}

// Compute derivatives for the 4 ODEs
void compute_derivatives(
    const BubbleState* state,
    const RPE_Params* params,
    BubbleDerivatives* derivs,
    CONVERGE_precision_t mdot
) {
    CONVERGE_precision_t R = state->R;
    CONVERGE_precision_t Rdot = state->Rdot;
    CONVERGE_precision_t T_drop = state->T_drop;
    CONVERGE_precision_t m_b = state->m_b;
    
    // Safety checks
    if (R < 1e-12) R = 1e-12;
    
    // Calculate bubble pressure
    CONVERGE_precision_t Vb = (4.0/3.0) * PI * R*R*R;
    CONVERGE_precision_t rho_v = safe_divide(m_b, Vb, 1e-6);
    if (rho_v < 1e-6) rho_v = 1e-6;
    
    CONVERGE_precision_t T_drop_safe = (T_drop > 1e-3) ? T_drop : 1e-3;
    CONVERGE_precision_t Pb = rho_v * params->R_spec * T_drop_safe;
    
    // 1. dR/dt = Rdot
    derivs->dRdt = Rdot;
    
    // 2. Rayleigh-Plesset equation: dRdot/dt
    CONVERGE_precision_t pressure_term = safe_divide(Pb - params->P_amb - 2.0*params->sigma/R - 4.0*params->mu_l*Rdot/R,
                                                       params->rho_l * R, 0.0);
    CONVERGE_precision_t inertial_term = -1.5 * Rdot * Rdot / R;
    CONVERGE_precision_t Rddot = pressure_term + inertial_term;
    
    // Apply acceleration limits
    if (Rddot > 1e10) Rddot = 1e10;
    if (Rddot < -1e10) Rddot = -1e10;
    
    derivs->dRdotdt = Rddot;
    
    // 3. Mass balance: dm_b/dt = mdot
    derivs->dmbdt = mdot;
    
    // 4. Energy balance: dT_drop/dt
    CONVERGE_precision_t Nu, Q_conv, mdot_local;
    compute_thermal_mass_transfer(R, Rdot, T_drop, m_b, params, &Nu, &Q_conv, &mdot_local);
    
    CONVERGE_precision_t Q_evap = params->L_v * mdot;
    CONVERGE_precision_t dTdt = safe_divide(Q_conv - Q_evap, params->m_drop * params->cp_l, 0.0);
    
    // Apply temperature rate limits
    if (dTdt > 1e6) dTdt = 1e6;
    if (dTdt < -1e6) dTdt = -1e6;
    
    derivs->dTdt = dTdt;
}

// Explicit Euler integration step
void euler_step(
    BubbleState* state,
    const BubbleDerivatives* derivs,
    CONVERGE_precision_t dt
) {
    state->R += derivs->dRdt * dt;
    state->Rdot += derivs->dRdotdt * dt;
    state->T_drop += derivs->dTdt * dt;
    state->m_b += derivs->dmbdt * dt;
    
    // Apply safety limits
    if (state->R < 1e-12) state->R = 1e-12;
    if (state->m_b < 0.0) state->m_b = 0.0;
    if (state->T_drop < 100.0) state->T_drop = 100.0;
    if (state->T_drop > 500.0) state->T_drop = 500.0;
}

// Main solver function - replaces Bubble_Velocity()
void RPE_euler_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
) {
    // Initialize parameters structure
    RPE_Params params;
    
    // Read liquid properties from parcel
    params.rho_l = old_parcel_cloud->density[p_idx];
    params.mu_l = old_parcel_cloud->viscosity[p_idx];
    params.sigma = old_parcel_cloud->surf_ten[p_idx];
    params.P_amb = P_amb;
    params.Ro = old_parcel_cloud->radius[p_idx];
    
    // Calculate droplet mass
    params.m_drop = (4.0/3.0) * PI * params.Ro * params.Ro * params.Ro * params.rho_l;
    
    // Get temperature-dependent properties from tables
    CONVERGE_precision_t Td = old_parcel_cloud->temp[p_idx];
    
    // Latent heat and specific heat (weighted by species mass fraction)
    CONVERGE_precision_t L_v_avg = 0.0;
    CONVERGE_precision_t cp_l_avg = 0.0;
    for (CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++) {
        L_v_avg += old_parcel_cloud->mfrac[p_idx * num_parcel_species + isp] * 
                   CONVERGE_table_lookup(hvap_table[isp], Td);
        cp_l_avg += old_parcel_cloud->mfrac[p_idx * num_parcel_species + isp] * 
                    CONVERGE_table_lookup(cp_table[isp], Td);
    }
    
    params.L_v = L_v_avg;
    params.cp_l = cp_l_avg;
    
    // Ammonia-specific properties (TODO: make these configurable)
    params.R_spec = R_SPEC_NH3;
    params.k_l = 0.5;  // Approximate thermal conductivity W/(m·K)
    params.max_Nu = 1000.0;  // Maximum Nusselt number
    
    // Initialize bubble state from parcel
    BubbleState state;
    state.R = old_parcel_cloud->r_bubble[p_idx];
    state.Rdot = old_parcel_cloud->v_bubble[p_idx];
    state.T_drop = Td;
    
    // Initialize bubble mass if not already set
    CONVERGE_precision_t P_sat;
    Saturation_PressureNH3(Td, &P_sat);
    CONVERGE_precision_t rho_b = bubble_densityNH3(P_sat, Td);
    CONVERGE_precision_t Vb = (4.0/3.0) * PI * state.R * state.R * state.R;
    state.m_b = rho_b * Vb;
    
    // Check for negative pressure difference
    if ((P_sat - P_amb) < 0.0) {
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        old_parcel_cloud->pbt[p_idx] = 0;
        old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
        return;
    }
    
    // Compute derivatives
    BubbleDerivatives derivs;
    CONVERGE_precision_t Nu, Q_conv, mdot;
    compute_thermal_mass_transfer(state.R, state.Rdot, state.T_drop, state.m_b, 
                                   &params, &Nu, &Q_conv, &mdot);
    compute_derivatives(&state, &params, &derivs, mdot);
    
    // Euler step
    euler_step(&state, &derivs, dt_sub);
    
    // Enforce constraint: R <= Ro
    if (state.R > params.Ro) {
        state.R = 0.95 * params.Ro;
        state.Rdot = 0.0;  // Stop growth at droplet edge
    }
    
    // Check for very small velocity
    if (state.Rdot < 1.0e-10) {
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        old_parcel_cloud->pbt[p_idx] = 0;
        old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
        return;
    }
    
    // Update parcel cloud with new state
    old_parcel_cloud->r_bubble[p_idx] = state.R;
    old_parcel_cloud->v_bubble[p_idx] = state.Rdot;
    old_parcel_cloud->temp[p_idx] = state.T_drop;
}
