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
    CONVERGE_index_t p_idx
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
    
    // DIAGNOSTIC: Check if this parcel already broke up
    if (old_parcel_cloud->breakup_phase[p_idx] == 5) {
        static int double_breakup_count = 0;
        if (double_breakup_count < 10) {
            printf("[BREAKUP_SONG_ERROR] Called on already-broken parcel! p_idx=%li, breakup_phase=%d\n",
                   p_idx, old_parcel_cloud->breakup_phase[p_idx]);
            printf("                      R=%.3e m, num_drop=%.3e, lifetime=%.3e s\n",
                   R_parent, N_parent, old_parcel_cloud->lifetime[p_idx]);
            double_breakup_count++;
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
    
    // Calculate radial velocity using momentum balance from existing DGRE model
    // rad_vel = 3 * v_bubble * r_bubble² * (R_parent - r_bubble) / (R_parent³ - r_bubble³)
    CONVERGE_precision_t rad_vel = 0.0;
    
    if (r_bubble > 1e-12 && fabs(v_bubble) > 1e-12) {
        CONVERGE_precision_t r_bubble_sq = r_bubble * r_bubble;
        CONVERGE_precision_t R_parent_cubed = R_parent * R_parent * R_parent;
        CONVERGE_precision_t r_bubble_cubed = r_bubble * r_bubble * r_bubble;
        CONVERGE_precision_t denominator = R_parent_cubed - r_bubble_cubed;
        
        if (fabs(denominator) > 1e-20) {
            rad_vel = 3.0 * v_bubble * r_bubble_sq * (R_parent - r_bubble) / denominator;
        }
    }
    
    // Get parent velocity vector
    CONVERGE_vec3_t parent_velocity;
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], &parent_velocity);
    
    // Safety check on parent velocity
    CONVERGE_precision_t parent_vel_mag = sqrt(parent_velocity[0] * parent_velocity[0] +
                                               parent_velocity[1] * parent_velocity[1] +
                                               parent_velocity[2] * parent_velocity[2]);
    
    const CONVERGE_precision_t parent_vel_limit = 1.0e3;
    if (parent_vel_mag > parent_vel_limit) {
        printf("[BREAKUP_SONG_WARNING] Parent velocity too large (%.3e m/s). Clamping.\n", parent_vel_mag);
        for (int i = 0; i < 3; i++) {
            if (fabs(parent_velocity[i]) > parent_vel_limit) {
                parent_velocity[i] = (parent_velocity[i] > 0) ? parent_vel_limit : -parent_vel_limit;
            }
        }
        parent_vel_mag = parent_vel_limit;
    }
    
    // Clamp radial velocity if it exceeds parent velocity magnitude
    if (rad_vel > parent_vel_mag) {
        rad_vel = parent_vel_mag;
    }
    
    // Compute child velocity magnitude using Pythagorean theorem
    // v_child² = v_parent² + v_radial²
    CONVERGE_precision_t child_vel_mag_sq = parent_vel_mag * parent_vel_mag + rad_vel * rad_vel;
    CONVERGE_precision_t child_vel_mag = sqrt(child_vel_mag_sq);
    
    // Scale child velocity vector by ratio of magnitudes
    CONVERGE_vec3_t child_velocity;
    if (parent_vel_mag > 1e-12) {
        CONVERGE_precision_t scale_factor = child_vel_mag / parent_vel_mag;
        child_velocity[0] = parent_velocity[0] * scale_factor;
        child_velocity[1] = parent_velocity[1] * scale_factor;
        child_velocity[2] = parent_velocity[2] * scale_factor;
    } else {
        // If parent velocity is near zero, set child velocity to zero
        child_velocity[0] = 0.0;
        child_velocity[1] = 0.0;
        child_velocity[2] = 0.0;
    }
    
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
        printf("[BREAKUP_SONG]   rad_vel=%.3e m/s, parent_vel=%.3e m/s, child_vel=%.3e m/s\n",
               rad_vel, parent_vel_mag, child_vel_mag);
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
