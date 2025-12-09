// Breakup_Song.c
// Song et al. breakup model - simplified approach
// Creates 2-5 child droplets with equal probability
// Modifies parent parcel in-place (no new parcel creation)
// Based on Song et al. paper breakup methodology

/*
 * breakup_phase states:
 *   0 = DISABLED  (parent, not eligible - subcooled, too small, etc.)
 *   1 = ELIGIBLE  (parent, superheated, ready to enter thermal breakup)
 *   2 = ACTIVE    (parent, growing bubble in sub-timestep loop)
 *   3 = RECOVERY  (parent, bubble collapsed, attempting recovery)
 *   4 = READY     (parent, bubble at threshold, ready to fragment)
 *   5 = COMPLETE  (child - result of breakup, any mechanism)
 */

#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Breakup_Song.h>
#include <PsatNH3.h>
#include <globals.h>

// Physical constants
#define PI 3.14159265358979323846

// Debug counter
static int song_breakup_logged = 0;
#define SONG_BREAKUP_LOG_MAX 10

// ============================================================================
// Song Breakup Function
// ============================================================================
// Randomly selects N_child ∈ {2, 3, 4, 5} with equal probability
// R_child = R_parent / cbrt(N_child_droplets)
// num_drop_child = num_drop_parent × N_child_droplets (volume balance)
// Updates parent parcel fields in-place (no new parcels created)
// ============================================================================
void Breakup_Song(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb
) {
    static int breakup_count = 0;
    
    // One-time logging to confirm Song breakup model is active
    if (breakup_count == 0) {
        printf("\n========================================\n");
        printf("[BREAKUP_SONG] Song breakup model active\n");
        printf("[BREAKUP_SONG] Creates 2-5 droplets with equal probability\n");
        printf("[BREAKUP_SONG] R_child = R_parent / cbrt(N_child)\n");
        printf("[BREAKUP_SONG] num_drop_child = num_drop_parent × N_child\n");
        printf("========================================\n\n");
    }
    
    // Read parent parcel properties
    CONVERGE_precision_t R_parent = old_parcel_cloud->radius[p_idx];
    CONVERGE_precision_t N_parent = old_parcel_cloud->num_drop[p_idx];
    CONVERGE_precision_t r_bubble = old_parcel_cloud->r_bubble[p_idx];
    CONVERGE_precision_t v_bubble = old_parcel_cloud->v_bubble[p_idx];
    
    // DIAGNOSTIC: Check if parcel is in correct state for breakup (must be READY = 4)
    if (old_parcel_cloud->breakup_phase[p_idx] != 4) {
        static int wrong_phase_count = 0;
        if (wrong_phase_count < 5) {
            printf("[BREAKUP_SONG_ERROR] Called on parcel NOT in READY state! p_idx=%li, breakup_phase=%d\n",
                   p_idx, old_parcel_cloud->breakup_phase[p_idx]);
            printf("                      R=%.3e m, num_drop=%.3e, lifetime=%.3e s\n",
                   R_parent, N_parent, old_parcel_cloud->lifetime[p_idx]);
            wrong_phase_count++;
        }
        // Still proceed with breakup but this indicates a logic error
    }
    
    // Safety checks
    if (R_parent < 1e-12 || N_parent < 1e-6) {
        printf("[BREAKUP_SONG_ERROR] Invalid parent: R=%.3e m, N=%.3e\n", R_parent, N_parent);
        return;
    }
    
    // CRITICAL CHECK: Prevent breakup of already-small parcels
    if (R_parent < 5e-6) {
        static int tiny_parent_count = 0;
        if (tiny_parent_count < 20) {
            printf("[BREAKUP_SONG_ABORT] Parent too small! R_parent=%.6e m (< 5 µm)\n", R_parent);
            printf("                      This indicates multiple breakups occurred!\n");
            printf("                      p_idx=%li, breakup_phase=%d, lifetime=%.3e s\n",
                   p_idx, old_parcel_cloud->breakup_phase[p_idx], old_parcel_cloud->lifetime[p_idx]);
            tiny_parent_count++;
        }
        // Still allow breakup for now, but this is a bug indicator
    }
    
    // Randomly select number of child droplets: 2, 3, 4, or 5 (equal probability)
    CONVERGE_precision_t rand_val = CONVERGE_random_precision();
    CONVERGE_index_t N_child_droplets;
    
    if (rand_val < 0.25) {
        N_child_droplets = 2;
    } else if (rand_val < 0.50) {
        N_child_droplets = 3;
    } else if (rand_val < 0.75) {
        N_child_droplets = 4;
    } else {
        N_child_droplets = 5;
    }
    
    // Calculate child droplet radius
    // When 1 droplet breaks into N pieces, each piece is smaller
    // R_child = R_parent / cbrt(N_child_droplets)
    CONVERGE_precision_t R_child = R_parent / cbrt((CONVERGE_precision_t)N_child_droplets);
    
    // Calculate child droplet number per parcel using volume balance
    // Conservation: num_drop_parent × R_parent³ = num_drop_child × R_child³
    // Therefore: num_drop_child = num_drop_parent × (R_parent/R_child)³ = num_drop_parent × N_child_droplets
    CONVERGE_precision_t N_child = N_parent * (CONVERGE_precision_t)N_child_droplets;
    
    // ========================================================================
    // ENERGY-BASED VELOCITY CALCULATION (Song et al. model)
    // ========================================================================
    
    // Get droplet temperature and calculate saturation pressure
    CONVERGE_precision_t T_drop = old_parcel_cloud->temp[p_idx];
    CONVERGE_precision_t P_sat;
    Saturation_PressureNH3(T_drop, &P_sat);
    
    // Get material properties from parcel cloud
    CONVERGE_precision_t sigma = old_parcel_cloud->surf_ten[p_idx];        // Surface tension
    CONVERGE_precision_t rho_l = old_parcel_cloud->density[p_idx];         // Liquid density
    CONVERGE_precision_t rho_g = old_parcel_cloud->gas_density[p_idx];     // Gas density
    
    // Calculate bubble volume
    CONVERGE_precision_t Vb = (4.0/3.0) * PI * r_bubble * r_bubble * r_bubble;
    
    // Calculate void fraction: ε = R_bubble³ / R_drop³
    CONVERGE_precision_t epsilon = (r_bubble * r_bubble * r_bubble) / (R_parent * R_parent * R_parent);
    if (epsilon > 1.0) epsilon = 1.0;
    if (epsilon < 0.0) epsilon = 0.0;
    
    // Calculate surface energy change
    // ΔE_surf = 4πσ(R² - (N1/N)*R1²)
    CONVERGE_precision_t surface_parent = 4.0 * PI * sigma * R_parent * R_parent;
    CONVERGE_precision_t surface_child = 4.0 * PI * sigma * R_child * R_child * (N_child / N_parent);
    CONVERGE_precision_t DeltaE_surf = surface_parent - surface_child;
    
    // Calculate pressure energy release
    // ΔE_p = (P_bubble - P_ambient) * V_bubble
    CONVERGE_precision_t DeltaE_p = (P_sat - P_amb) * Vb;
    
    // Calculate mixture density factor (from Song paper: 0.45*(ε*ρ_g + (1-ε)*ρ_l))
    CONVERGE_precision_t rho_m = 0.45 * (epsilon * rho_g + (1.0 - epsilon) * rho_l);
    
    // Safety check on mixture density
    if (rho_m < 1e-6) {
        rho_m = rho_l * 0.45;  // Fallback to liquid density
    }
    
    // Calculate velocity increment magnitude
    // ΔU = sqrt(2*(ΔE_surf + ΔE_p) / (ρ_m * V_b))
    CONVERGE_precision_t energy_sum = DeltaE_surf + DeltaE_p;
    CONVERGE_precision_t DeltaU = 0.0;
    
    if (energy_sum > 0.0 && Vb > 1e-30 && rho_m > 1e-6) {
        DeltaU = sqrt(2.0 * energy_sum / (rho_m * Vb));
    }
    
    // Safety clamp on velocity increment (max 500 m/s)
    const CONVERGE_precision_t max_delta_v = 500.0;
    if (DeltaU > max_delta_v) {
        DeltaU = max_delta_v;
    }
    
    // Get parent velocity vector
    CONVERGE_vec3_t parent_velocity;
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], &parent_velocity);
    
    // Calculate parent velocity magnitude
    CONVERGE_precision_t ux = parent_velocity[0];
    CONVERGE_precision_t uy = parent_velocity[1];
    CONVERGE_precision_t uz = parent_velocity[2];
    CONVERGE_precision_t ud = sqrt(ux*ux + uy*uy + uz*uz);
    
    // Calculate directional velocity increments (along parent velocity direction)
    CONVERGE_precision_t DeltaUx = 0.0;
    CONVERGE_precision_t DeltaUy = 0.0;
    CONVERGE_precision_t DeltaUz = 0.0;
    
    if (ud > 1e-12) {
        // Direction is along parent velocity
        DeltaUx = (ux / ud) * DeltaU;
        DeltaUy = (uy / ud) * DeltaU;
        DeltaUz = (uz / ud) * DeltaU;
    }
    
    // Calculate child velocity by adding increment to parent velocity
    CONVERGE_vec3_t child_velocity;
    child_velocity[0] = ux + DeltaUx;
    child_velocity[1] = uy + DeltaUy;
    child_velocity[2] = uz + DeltaUz;
    
    // Calculate child velocity magnitude for diagnostics
    CONVERGE_precision_t child_vel_mag = sqrt(child_velocity[0]*child_velocity[0] + 
                                              child_velocity[1]*child_velocity[1] + 
                                              child_velocity[2]*child_velocity[2]);
    
    // Diagnostic logging for ALL breakups (increased from 10 to 50)
    static int total_breakup_count = 0;
    total_breakup_count++;
    
    if (song_breakup_logged < 50) {  // Increased logging
        printf("\n[BREAKUP_SONG #%d] ==========================================\n", total_breakup_count);
        printf("[BREAKUP_SONG] p_idx=%li: N_child_droplets=%li, lifetime=%.3e s\n", 
               p_idx, N_child_droplets, old_parcel_cloud->lifetime[p_idx]);
        printf("[BREAKUP_SONG]   BEFORE: R=%.6e m, num_drop=%.6e\n", R_parent, N_parent);
        printf("[BREAKUP_SONG]   AFTER:  R=%.6e m, num_drop=%.6e\n", R_child, N_child);
        printf("[BREAKUP_SONG]   r_bubble=%.3e m, v_bubble=%.3e m/s\n", r_bubble, v_bubble);
        printf("[BREAKUP_SONG]   ENERGY: DeltaE_surf=%.3e J, DeltaE_p=%.3e J\n", DeltaE_surf, DeltaE_p);
        printf("[BREAKUP_SONG]   MIXTURE: epsilon=%.4f, rho_m=%.3e kg/m³\n", epsilon, rho_m);
        printf("[BREAKUP_SONG]   PRESSURE: P_sat=%.3e Pa, P_amb=%.3e Pa, dP=%.3e Pa\n", P_sat, P_amb, P_sat - P_amb);
        printf("[BREAKUP_SONG]   VELOCITY: parent=%.3e m/s, DeltaU=%.3e m/s, child=%.3e m/s\n",
               ud, DeltaU, child_vel_mag);
        printf("[BREAKUP_SONG]   FLAGS BEFORE: breakup_phase=%d\n",
               old_parcel_cloud->breakup_phase[p_idx]);
        
        // Verify volume conservation
        // Volume balance: num_drop_parent × R_parent³ = num_drop_child × R_child³
        CONVERGE_precision_t vol_parent = N_parent * R_parent * R_parent * R_parent;
        CONVERGE_precision_t vol_child = N_child * R_child * R_child * R_child;
        CONVERGE_precision_t vol_error = fabs(vol_child - vol_parent) / (vol_parent + 1e-30);
        printf("[BREAKUP_SONG]   Volume conservation error: %.2e%%\n", vol_error * 100.0);
        
        // Calculate reduction ratios
        printf("[BREAKUP_SONG]   Radius ratio: %.4f (child/parent)\n", R_child/R_parent);
        printf("[BREAKUP_SONG]   num_drop ratio: %.4f (child/parent)\n", N_child/N_parent);
        
        song_breakup_logged++;
    }
    
    // Update parent parcel in-place to become a child parcel
    // NOTE: This represents ONE of the N_child_droplets created
    // The other (N_child_droplets - 1) droplets are represented implicitly by reducing num_drop
    
    old_parcel_cloud->radius[p_idx] = R_child;
    old_parcel_cloud->num_drop[p_idx] = N_child;
    
    // Update velocity
    old_parcel_cloud->uu[p_idx][0] = child_velocity[0];
    old_parcel_cloud->uu[p_idx][1] = child_velocity[1];
    old_parcel_cloud->uu[p_idx][2] = child_velocity[2];
    
    // Mark as child parcel (phase 5) and reset bubble
    old_parcel_cloud->breakup_phase[p_idx] = 5;  // COMPLETE (child)
    old_parcel_cloud->film_flag[p_idx] = 5;  // Hijack: mirror breakup_phase
    old_parcel_cloud->r_bubble[p_idx] = 0.0;
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    
    // DIAGNOSTIC: Confirm phase was set
    if (song_breakup_logged <= 50) {
        printf("[BREAKUP_SONG]   FLAGS AFTER:  breakup_phase=%d\n",
               old_parcel_cloud->breakup_phase[p_idx]);
        printf("[BREAKUP_SONG] ==========================================\n\n");
    }
    
    // Reset bubble properties (bubble is destroyed during breakup)
    old_parcel_cloud->r_bubble[p_idx] = 0.0;
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    old_parcel_cloud->r_bubble_0[p_idx] = 0.0;
    
    // Update mass-related fields
    // m0 = (4/3) * π * R³ * num_drop (liquid mass only)
    old_parcel_cloud->m0[p_idx] = (4.0/3.0) * PI * R_child * R_child * R_child * N_child;
    
    // DIAGNOSTIC: Final state summary
    if (song_breakup_logged <= 50) {
        CONVERGE_precision_t m0_parent = (4.0/3.0) * PI * R_parent * R_parent * R_parent * N_parent;
        CONVERGE_precision_t m0_child = old_parcel_cloud->m0[p_idx];
        printf("[BREAKUP_SONG] FINAL STATE: R=%.6e m, num_drop=%.6e, m0=%.6e\n", 
               old_parcel_cloud->radius[p_idx], old_parcel_cloud->num_drop[p_idx], m0_child);
        printf("[BREAKUP_SONG] m0 conservation: parent=%.6e, child=%.6e, error=%.2e%%\n\n",
               m0_parent, m0_child, fabs(m0_child - m0_parent)/m0_parent * 100.0);
    }
    
    // Optional: Update film_flag for diagnostic tracking (if not using film model)
    // old_parcel_cloud->film_flag[p_idx] = 1;  // Uncomment if needed for diagnostics
    
    breakup_count++;
}
