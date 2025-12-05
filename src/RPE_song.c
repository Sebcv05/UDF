//RPE_song.c
// Song et al. Rayleigh-Plesset Equation solver (isothermal model)
// Based on standalone test_song_rpe.c, adapted for CONVERGE UDF
//
// IMPORTANT: This solver handles bubble GROWTH only (void fraction criterion).
// The DGRE model (called separately in spray_drop_distort_NH3.c) handles bubble
// geometry evolution and determines final breakup based on its own criteria.
//
// Void fraction = 0.55 triggers reset_parcel_to_child, after which DGRE takes over.

#include "lagrangian/env.h"
#include <RPE_song.h>
#include <PsatNH3.h>
#include <parcel_reset.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>

// Physical constants
#define PI 3.14159265358979323846
#define R_SPEC_NH3 488.2  // J/(kg·K) - Specific gas constant for ammonia

// Song model parameters
#define P_R0_SONG 1.0e6   // Residual gas pressure (Pa)
#define KAPPA_SONG 0.0    // Surface viscosity (negligible)
#define VOID_TARGET 0.55  // Termination criterion: bubble reaches breakup threshold
                          // After this, DGRE model determines final breakup geometry

// Debug logging
static int song_debug_count = 0;
#define SONG_DEBUG_MAX 10

// ============================================================================
// Helper Function: Compute void fraction
// ============================================================================
// Correct formula: ε = V_bubble / (V_bubble + V_liquid)
//                    = R_bubble³ / (R_bubble³ + R_liquid³)
// 
// For a droplet: V_total = (4/3)π·R_drop³
//                V_bubble = (4/3)π·R_bubble³
//                V_liquid = V_total - V_bubble
//
// Therefore: R_liquid³ = R_drop³ - R_bubble³
//            ε = R_bubble³ / R_drop³
CONVERGE_precision_t song_compute_void_fraction(
    CONVERGE_precision_t R_bubble,
    CONVERGE_precision_t R_drop
) {
    // ε = R_bubble³ / R_drop³
    CONVERGE_precision_t R_bubble_cubed = R_bubble * R_bubble * R_bubble;
    CONVERGE_precision_t R_drop_cubed = R_drop * R_drop * R_drop;
    
    if (R_drop_cubed < 1e-30) return 0.0;
    
    CONVERGE_precision_t epsilon = R_bubble_cubed / R_drop_cubed;
    
    // Clamp to [0, 1]
    if (epsilon < 0.0) epsilon = 0.0;
    if (epsilon > 1.0) epsilon = 1.0;
    
    return epsilon;
}

// ============================================================================
// Helper Function: Compute mixture density
// ============================================================================
CONVERGE_precision_t song_compute_mixture_density(
    CONVERGE_precision_t epsilon,
    CONVERGE_precision_t rho_v,
    CONVERGE_precision_t rho_l
) {
    // ρ_m = ε·ρ_v + (1-ε)·ρ_l
    return epsilon * rho_v + (1.0 - epsilon) * rho_l;
}

// ============================================================================
// Helper Function: Compute Song RPE acceleration
// ============================================================================
// Song RPE (Equation 4):
// (3/2)·ρ_m·Ṙ² + ρ_m·R·R̈ = P_sat - P_∞ + (2σ/R₀+P_r0)·(R₀/R)³ - 2σ/R - 4μ·Ṙ/R - 4κ·Ṙ/R²
//
// Solving for R̈:
// R̈ = [P_sat - P_∞ + (2σ/R₀+P_r0)·(R₀/R)³ - 2σ/R - 4μ·Ṙ/R - 4κ·Ṙ/R²] / (ρ_m·R) - (3/2)·Ṙ²/R
// ============================================================================
CONVERGE_precision_t song_compute_acceleration(
    CONVERGE_precision_t R,
    CONVERGE_precision_t Rdot,
    CONVERGE_precision_t R0,
    CONVERGE_precision_t rho_m,
    CONVERGE_precision_t P_sat,
    const SongParams* params
) {
    // Prevent division by zero
    if (R < 1e-12) R = 1e-12;
    if (rho_m < 1e-6) rho_m = 1e-6;
    if (R0 < 1e-12) R0 = 1e-12;
    
    // Initial pressure term: (2σ/R0 + P_r0) * (R0/R)³
    CONVERGE_precision_t P_init_coeff = 2.0 * params->sigma / R0 + params->P_r0;
    CONVERGE_precision_t R_ratio = R0 / R;
    CONVERGE_precision_t R_ratio_cubed = R_ratio * R_ratio * R_ratio;
    CONVERGE_precision_t P_init = P_init_coeff * R_ratio_cubed;
    
    // Pressure terms: P_sat - P_∞ + P_init - 2σ/R
    CONVERGE_precision_t pressure_term = P_sat - params->P_amb + P_init - 2.0 * params->sigma / R;
    
    // Viscous dissipation: -4μ·Ṙ/R - 4κ·Ṙ/R²
    CONVERGE_precision_t viscous_term = -4.0 * params->mu_l * Rdot / R 
                                       - 4.0 * params->kappa * Rdot / (R * R);
    
    // Combined numerator
    CONVERGE_precision_t numerator = pressure_term + viscous_term;
    
    // Inertial term: -(3/2)·Ṙ²/R
    CONVERGE_precision_t inertial_term = -1.5 * Rdot * Rdot / R;
    
    // R̈ = numerator / (ρ_m·R) + inertial
    CONVERGE_precision_t Rddot = numerator / (rho_m * R) + inertial_term;
    
    return Rddot;
}

// ============================================================================
// Main Song RPE Solver
// ============================================================================
void RPE_song_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
) {
    // One-time logging to confirm Song model is active
    static int song_model_logged = 0;
    if (!song_model_logged) {
        printf("\n========================================\n");
        printf("[RPE_MODEL_SONG] Song isothermal RPE model active\n");
        printf("========================================\n\n");
        song_model_logged = 1;
    }
    
    // Initialize parameters from parcel properties
    SongParams params;
    params.rho_l = old_parcel_cloud->density[p_idx];
    params.mu_l = old_parcel_cloud->viscosity[p_idx];
    params.sigma = old_parcel_cloud->surf_ten[p_idx];
    params.R_spec = R_SPEC_NH3;
    params.P_amb = P_amb;
    params.P_r0 = P_R0_SONG;
    params.kappa = KAPPA_SONG;
    
    // Read current state from parcel
    CONVERGE_precision_t R = old_parcel_cloud->r_bubble[p_idx];
    CONVERGE_precision_t Rdot = old_parcel_cloud->v_bubble[p_idx];
    CONVERGE_precision_t R0 = old_parcel_cloud->r_bubble_0[p_idx];
    CONVERGE_precision_t Ro = old_parcel_cloud->r_drop_0[p_idx];
    CONVERGE_precision_t R_drop_current = old_parcel_cloud->radius[p_idx];  // Current droplet radius
    CONVERGE_precision_t T_drop = old_parcel_cloud->temp[p_idx];
    
    // Safety checks on initial values
    if (Ro < 1e-9) {
        printf("[SONG_ERROR] Droplet radius too small: Ro=%.3e m\n", Ro);
        // Don't reset - just return, main loop will handle
        return;
    }
    
    if (R0 < 1e-12) {
        printf("[SONG_ERROR] Initial bubble radius too small: R0=%.3e m\n", R0);
        // Don't reset - just return, main loop will handle
        return;
    }
    
    // Calculate saturation pressure (isothermal - temperature is constant)
    CONVERGE_precision_t P_sat;
    Saturation_PressureNH3(T_drop, &P_sat);
    
    // Check superheat condition
    if (P_sat <= P_amb) {
        if (song_debug_count < SONG_DEBUG_MAX) {
            printf("[SONG_ABORT] Not superheated: P_sat=%.2e Pa, P_amb=%.2e Pa\n", P_sat, P_amb);
            song_debug_count++;
        }
        // Don't reset - just return, main loop will handle
        return;
    }
    
    // Compute void fraction using CURRENT droplet radius
    CONVERGE_precision_t epsilon = song_compute_void_fraction(R, R_drop_current);
    
    // NOTE: Void fraction check moved to main loop in spray_drop_distort_NH3.c
    // Song solver just grows bubble, main loop checks epsilon and triggers breakup
    
    // Compute vapor density (ideal gas law)
    CONVERGE_precision_t rho_v = P_sat / (params.R_spec * T_drop);
    if (rho_v < 1e-6) rho_v = 1e-6;
    
    // Compute mixture density
    CONVERGE_precision_t rho_m = song_compute_mixture_density(epsilon, rho_v, params.rho_l);
    
    // Debug logging for first few calls
    if (song_debug_count < SONG_DEBUG_MAX) {
        printf("[SONG_STEP] p_idx=%li, R=%.3e m, Rdot=%.3e m/s, ε=%.4f, ρ_m=%.1f kg/m³, T=%.2f K\n",
               p_idx, R, Rdot, epsilon, rho_m, T_drop);
        printf("[SONG_STEP]   R0=%.3e m, P_sat=%.2e Pa, P_amb=%.2e Pa, ΔP=%.2e Pa\n",
               R0, P_sat, P_amb, P_sat - P_amb);
        song_debug_count++;
    }
    
    // Compute acceleration
    CONVERGE_precision_t Rddot = song_compute_acceleration(R, Rdot, R0, rho_m, P_sat, &params);
    
    // Explicit Euler integration (single timestep)
    Rdot += Rddot * dt_sub;
    R += Rdot * dt_sub;
    
    // Safety checks
    if (R < 1e-12) R = 1e-12;
    
    // Check for collapse (negative velocity)
    if (Rdot < 0.0) {
        static int collapse_detail_count = 0;
        if (collapse_detail_count < 5) {
            printf("[SONG_COLLAPSE] Bubble collapsing: Rdot=%.3e m/s at p_idx=%li\n", Rdot, p_idx);
            printf("[SONG_COLLAPSE]   R=%.3e m, R0=%.3e m, Rddot=%.3e m/s²\n", R, R0, Rddot);
            printf("[SONG_COLLAPSE]   P_sat=%.2e Pa, P_amb=%.2e Pa, ΔP=%.2e Pa\n", P_sat, P_amb, P_sat - P_amb);
            printf("[SONG_COLLAPSE]   T_drop=%.2f K, rho_m=%.1f kg/m³, ε=%.4f\n", T_drop, rho_m, epsilon);
            
            // Calculate individual pressure terms for diagnosis
            CONVERGE_precision_t P_init_coeff = 2.0 * params.sigma / R0 + params.P_r0;
            CONVERGE_precision_t R_ratio_cubed = pow(R0/R, 3.0);
            CONVERGE_precision_t P_init = P_init_coeff * R_ratio_cubed;
            CONVERGE_precision_t P_laplace = 2.0 * params.sigma / R;
            printf("[SONG_COLLAPSE]   P_init=%.2e Pa, P_laplace=%.2e Pa\n", P_init, P_laplace);
            collapse_detail_count++;
        }
        // Don't reset - just return, main loop will handle
        return;
    }
    
    // Check if bubble reached droplet edge
    if (R > 0.999 * Ro) {
        R = 0.999 * Ro;
        Rdot = 0.0;
        printf("[SONG_EDGE] Bubble capped at droplet edge: R=%.3e m, Ro=%.3e m\n", R, Ro);
    }
    
    // Update parcel state
    // IMPORTANT: Temperature is NOT updated (isothermal model)
    old_parcel_cloud->r_bubble[p_idx] = R;
    old_parcel_cloud->v_bubble[p_idx] = Rdot;
    // old_parcel_cloud->temp[p_idx] remains unchanged
}
