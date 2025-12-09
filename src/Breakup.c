// Breakup.c

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
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <Breakup.h>
#include <spray_break.h>
#include <PsatNH3.h>
#include <globals.h>
#include <Vb.h>

// Global variables
// Global velocity variables
CONVERGE_precision_t user_child_velocity_x =0.0;
CONVERGE_precision_t user_child_velocity_y = 0.0;
CONVERGE_precision_t user_child_velocity_z = 0.0;
CONVERGE_precision_t breakup_velocity_scale = 1.0;
CONVERGE_precision_t breakup_radius_scale = 1.0;
CONVERGE_precision_t kb_threshold = 1.0;
CONVERGE_index_t num_child_parcels = 12;
static int breakup_scale_logged = 0;

// Langmuir-Knudsen evaporation model parameters (set by read_input.c)
CONVERGE_index_t lk_correction_flag = 0;
CONVERGE_index_t lk_diagnostic_flag = 0;
CONVERGE_precision_t lk_chi_neq_min = 0.0;
CONVERGE_precision_t lk_chi_neq_max = 0.9999;

// Song RPE model selection (default: thermal model)
CONVERGE_index_t use_song_rpe = 0;

// Profiling accumulators
static double prof_calcs = 0.0;
static double prof_loop = 0.0;
static double prof_child_parcel = 0.0;
static double prof_property_copy = 0.0;
static double prof_velocity_calc = 0.0;
static double prof_insert_cloud = 0.0;
static int last_cycle = -1;

// Rosin-Rammler global parameters (initialized once)
static RR_Params rr_params = {0};

// ============================================================================
// Gamma function approximation (Lanczos method)
// Only called during initialization, not per-breakup
// ============================================================================
static double tgamma_lanczos(double z) {
    static const double g = 7.0;
    static const double coef[9] = {
        0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
    };
    
    if (z < 0.5) {
        return M_PI / (sin(M_PI * z) * tgamma_lanczos(1.0 - z));
    }
    
    z -= 1.0;
    double x = coef[0];
    for (int i = 1; i < 9; i++) {
        x += coef[i] / (z + i);
    }
    
    double t = z + g + 0.5;
    return sqrt(2.0 * M_PI) * pow(t, z + 0.5) * exp(-t) * x;
}

// ============================================================================
// Initialize Rosin-Rammler parameters (called once at startup)
// ============================================================================
void init_RR_distribution(double n_shape) {
    if (rr_params.initialized) {
        return;
    }
    
    rr_params.n_RR = n_shape;
    
    // Pre-compute gamma ratio: tgamma(1+2/n) / tgamma(1+3/n)
    double gamma_2n = tgamma_lanczos(1.0 + 2.0 / n_shape);
    double gamma_3n = tgamma_lanczos(1.0 + 3.0 / n_shape);
    rr_params.gamma_ratio = gamma_2n / gamma_3n;
    
    rr_params.initialized = 1;
    
    int rank;
    CONVERGE_mpi_comm_rank(&rank);
    if (rank == 0) {
        printf("[RR_INIT] Rosin-Rammler distribution initialized:\n");
        printf("[RR_INIT]   n_RR = %.3f\n", rr_params.n_RR);
        printf("[RR_INIT]   tgamma(1+2/n) = %.6e\n", gamma_2n);
        printf("[RR_INIT]   tgamma(1+3/n) = %.6e\n", gamma_3n);
        printf("[RR_INIT]   gamma_ratio = %.6e\n", rr_params.gamma_ratio);
    }
}

// ============================================================================
// FALLBACK: Uniform distribution (used if RR sampling fails)
// ============================================================================
void fallback_uniform_children(
    CONVERGE_precision_t parent_radius,
    CONVERGE_precision_t parent_num_drop,
    CONVERGE_precision_t R32_target,
    int N,
    CONVERGE_precision_t* child_radii,
    CONVERGE_precision_t* child_num_drop
) {
    // All children get same radius (R32_target) and same num_drop
    CONVERGE_precision_t parent_volume = CONVERGE_cube(parent_radius);
    CONVERGE_precision_t child_volume = CONVERGE_cube(R32_target);
    CONVERGE_precision_t uniform_num_drop = parent_num_drop * parent_volume / (N * child_volume);
    
    for (int i = 0; i < N; i++) {
        child_radii[i] = R32_target;
        child_num_drop[i] = uniform_num_drop;
    }
    
    static int fallback_warn_count = 0;
    if (fallback_warn_count < 3) {
        printf("[RR_FALLBACK] Using uniform distribution: all R=%.3e m, num_drop=%.3e\n",
               R32_target, uniform_num_drop);
        fallback_warn_count++;
    }
}

// ============================================================================
// ROSIN-RAMMLER CHILD SAMPLING FUNCTION
// Returns: 0 on success, -1 on failure (fallback to uniform distribution)
// ============================================================================
int sample_RR_children(
    CONVERGE_precision_t parent_radius,
    CONVERGE_precision_t parent_num_drop,
    CONVERGE_precision_t R32_target,
    int N,
    double n_RR,
    double gamma_ratio,
    CONVERGE_precision_t* child_radii,
    CONVERGE_precision_t* child_num_drop
) {
    const double PI = 3.14159265358979323846;
    const double U_MIN = 1.0e-12;
    const double U_MAX = 1.0 - 1.0e-12;
    const double SCALE_MIN = 0.1;
    const double SCALE_MAX = 10.0;
    
    static int error_count = 0;
    const int MAX_ERRORS = 5;
    
    // Input validation
    if (parent_radius <= 0.0 || parent_num_drop <= 0.0 || R32_target <= 0.0 || N <= 0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Invalid inputs: R_p=%.3e, Nd_p=%.3e, R32=%.3e, N=%d\n",
                   parent_radius, parent_num_drop, R32_target, N);
            error_count++;
        }
        return -1;
    }
    
    if (n_RR <= 0.0 || gamma_ratio <= 0.0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Invalid RR parameters: n_RR=%.3f, gamma_ratio=%.6f\n",
                   n_RR, gamma_ratio);
            error_count++;
        }
        return -1;
    }
    
    CONVERGE_precision_t D32_target = 2.0 * R32_target;
    CONVERGE_precision_t X_RR = D32_target * gamma_ratio;
    double inv_n_RR = 1.0 / n_RR;
    
    CONVERGE_precision_t child_diameters[N];
    CONVERGE_precision_t total_sampled_volume = 0.0;
    
    for (int i = 0; i < N; i++) {
        CONVERGE_precision_t u = CONVERGE_random_precision();
        if (u < U_MIN) u = U_MIN;
        if (u > U_MAX) u = U_MAX;
        
        double log_arg = 1.0 - u;
        if (log_arg <= 0.0) {
            if (error_count < MAX_ERRORS) {
                printf("[RR_ERROR] Invalid log argument: u=%.12e\n", u);
                error_count++;
            }
            return -1;
        }
        
        child_diameters[i] = X_RR * pow(-log(log_arg), inv_n_RR);
        child_radii[i] = child_diameters[i] / 2.0;
        total_sampled_volume += CONVERGE_cube(child_radii[i]);
    }
    
    if (total_sampled_volume <= 0.0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Zero sampled volume: %.3e\n", total_sampled_volume);
            error_count++;
        }
        return -1;
    }
    
    CONVERGE_precision_t D2_sum = 0.0, D3_sum = 0.0;
    for (int i = 0; i < N; i++) {
        CONVERGE_precision_t D = child_diameters[i];
        D2_sum += D * D;
        D3_sum += D * D * D;
    }
    
    if (D2_sum <= 0.0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Zero D2_sum\n");
            error_count++;
        }
        return -1;
    }
    
    CONVERGE_precision_t D32_sample = D3_sum / D2_sum;
    if (D32_sample <= 0.0 || !isfinite(D32_sample)) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Invalid D32_sample: %.3e\n", D32_sample);
            error_count++;
        }
        return -1;
    }
    
    CONVERGE_precision_t scale_correction = D32_target / D32_sample;
    if (scale_correction < SCALE_MIN || scale_correction > SCALE_MAX) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_WARNING] Extreme scale: %.4f (clamping)\n", scale_correction);
            error_count++;
        }
        if (scale_correction < SCALE_MIN) scale_correction = SCALE_MIN;
        if (scale_correction > SCALE_MAX) scale_correction = SCALE_MAX;
    }
    
    total_sampled_volume = 0.0;
    for (int i = 0; i < N; i++) {
        child_diameters[i] *= scale_correction;
        child_radii[i] = child_diameters[i] / 2.0;
        total_sampled_volume += CONVERGE_cube(child_radii[i]);
    }
    
    if (total_sampled_volume <= 0.0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Zero volume after scaling\n");
            error_count++;
        }
        return -1;
    }
    
    CONVERGE_precision_t parent_volume = CONVERGE_cube(parent_radius);
    CONVERGE_precision_t base_num_drop = parent_num_drop * parent_volume / total_sampled_volume;
    
    if (base_num_drop <= 0.0 || !isfinite(base_num_drop)) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Invalid base_num_drop: %.3e\n", base_num_drop);
            error_count++;
        }
        return -1;
    }
    
    for (int i = 0; i < N; i++) {
        child_num_drop[i] = base_num_drop;
    }
    
    static int rr_diag_count = 0;
    if (rr_diag_count < 3) {
        // printf("[RR_SAMPLE] Parent: R=%.3e m, num_drop=%.3e\n", parent_radius, parent_num_drop);
        // printf("[RR_SAMPLE] Target R32=%.3e m, X_RR=%.3e m, n_RR=%.2f\n", R32_target, X_RR, n_RR);
        // printf("[RR_SAMPLE] D32_sample=%.3e (ratio=%.4f), scale=%.4f\n", 
        //        D32_sample, D32_sample/D32_target, scale_correction);
        // printf("[RR_SAMPLE] base_num_drop=%.3e (same for all)\n", base_num_drop);
        
        CONVERGE_precision_t total_child_volume = 0.0;
        for (int i = 0; i < N; i++) {
            total_child_volume += child_num_drop[i] * CONVERGE_cube(child_radii[i]);
        }
        CONVERGE_precision_t volume_error = fabs(total_child_volume - parent_num_drop * parent_volume) / 
                                           (parent_num_drop * parent_volume);
        // printf("[RR_SAMPLE] Volume conservation: error=%.2e%% (num_drop*R^3, no 4/3*pi)\n", volume_error * 100.0);
        
        if (volume_error > 1.0e-6) {
            printf("[RR_WARNING] Volume error > 1 ppm\n");
        }
        rr_diag_count++;
    }
    
    return 0;
}

// Function to print profiling information


void Breakup(struct ParcelCloud *old_parcel_cloud, CONVERGE_index_t p_idx, CONVERGE_cloud_t cloud)
{
    if(!breakup_scale_logged)
    {
        breakup_scale_logged = 1;
    }

    // UDF-level check to prevent creating child parcels from a parent with non-physical properties
    if (isnan(old_parcel_cloud->temp[p_idx]) || isinf(old_parcel_cloud->temp[p_idx]) || old_parcel_cloud->temp[p_idx] < 100.0) {
        CONVERGE_logger_warn("Breakup.c: Parent parcel %d (cloud %d) has invalid temperature (%.2f) at ncyc %ld. Skipping breakup for this parcel.",
            p_idx, old_parcel_cloud->cloud_index[p_idx], CONVERGE_ncyc(), old_parcel_cloud->temp[p_idx]);
        return;
    }
    if (isnan(old_parcel_cloud->uu[p_idx][0]) || isinf(old_parcel_cloud->uu[p_idx][0]) ||
        isnan(old_parcel_cloud->uu[p_idx][1]) || isinf(old_parcel_cloud->uu[p_idx][1]) ||
        isnan(old_parcel_cloud->uu[p_idx][2]) || isinf(old_parcel_cloud->uu[p_idx][2])) {
        CONVERGE_logger_warn("Breakup.c: Parent parcel %d (cloud %d) has invalid velocity at ncyc %ld. Skipping breakup for this parcel.",
            p_idx, old_parcel_cloud->cloud_index[p_idx], CONVERGE_ncyc());
        return;
    }


    //Timing vars
    CONVERGE_precision_t t0 = CONVERGE_mpi_wtime();
    // Section 1: Initialization and validation
    
    // Check if cloud and parcel cloud are valid
    if (!cloud || !old_parcel_cloud) {
        printf("\nBreakup.c: Invalid cloud or parcel cloud pointer\n");
        CONVERGE_mpi_abort();
    }
    CONVERGE_vec3_t user_child_velocity[MAX_NUM_CHILDREN];

    // Get cloud size and verify parcel index
    CONVERGE_index_t cloud_size = CONVERGE_cloud_size(cloud);
    // printf("\nBreakup.c: Cloud size = %d, p_idx = %d\n", cloud_size, p_idx);
    if (p_idx >= cloud_size) {
        printf("\nBreakup.c: Invalid parcel index %ld (cloud size = %ld)\n", p_idx, cloud_size);
        CONVERGE_mpi_abort();
    }

    // Verify all required fields are loaded
    if (!old_parcel_cloud->child_uu) {
        printf("\nBreakup.c: child_uu field not loaded\n");
        CONVERGE_mpi_abort();
    }
    if (!old_parcel_cloud->uu) {
        printf("\nBreakup.c: uu field not loaded\n");
        CONVERGE_mpi_abort();
    }
    if (!old_parcel_cloud->radius) {
        printf("\nBreakup.c: radius field not loaded\n");
        CONVERGE_mpi_abort();
    }
    if (!old_parcel_cloud->r_bubble) {
        printf("\nBreakup.c: r_bubble field not loaded\n");
        CONVERGE_mpi_abort();
    }
    if (!old_parcel_cloud->v_bubble) {
        printf("\nBreakup.c: v_bubble field not loaded\n");
        CONVERGE_mpi_abort();
    }

    // Verify parent velocity exists
    if (!old_parcel_cloud->uu[p_idx]) {
        printf("\nBreakup.c: Parent velocity at p_idx %ld is NULL\n", p_idx);
        CONVERGE_mpi_abort();
    }
    // printf("\nBreakup.c: Parent velocity = %e %e %e\n", 
        //    old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2]);

    // Verify child_uu exists for this parcel
    if (!old_parcel_cloud->child_uu[p_idx]) {
        printf("\nBreakup.c: child_uu at p_idx %ld is NULL\n", p_idx);
        CONVERGE_mpi_abort();
    }
    if(old_parcel_cloud->breakup_phase[p_idx] == 5){
    printf("\n ERROR, breakup routine being triggered on child parcel (breakup_phase = 5), phase = %i",old_parcel_cloud->breakup_phase[p_idx]);
    CONVERGE_mpi_abort();
}
 CONVERGE_precision_t   old_r = old_parcel_cloud->radius[p_idx];
 CONVERGE_precision_t  old_nd = old_parcel_cloud->num_drop[p_idx];
 CONVERGE_precision_t  old_r_bubble = old_parcel_cloud->r_bubble[p_idx];
 
    // DIAGNOSTIC: Print parameters EVERY TIME breakup runs
    printf("[PARAM_CHECK] breakup_velocity_scale=%.3f radius_scale=%.3f kb_threshold=%.3f num_children=%ld\n",
           breakup_velocity_scale, breakup_radius_scale, kb_threshold, num_child_parcels);
    fflush(stdout);
 
    // DIAGNOSTIC: Check for large parcels entering breakup
    static FILE* breakup_log = NULL;
    if (!breakup_log) {
        breakup_log = fopen("breakup_debug.csv", "w");
        if (breakup_log) {
            fprintf(breakup_log, "time,ncyc,p_idx,event,R_drop_before,R_bubble_before,R_drop_after,R_bubble_after,num_drop_before,num_drop_after,temp,film_flag,mean_child_radius\n");
        }
    }
    
    // CSV logger for breakup event data
    static FILE* breakup_events_log = NULL;
    if (!breakup_events_log) {
        breakup_events_log = fopen("breakup_events.csv", "w");
        if (breakup_events_log) {
            fprintf(breakup_events_log, "time,ncyc,p_idx,parent_radius,r_bubble,v_bubble,child_radius\n");
        }
    }
    
    // DON'T log entry yet - wait until we know child radii
    
    // printf("\nbreakup count %i",breakup_counter);
    if (old_parcel_cloud->radius[p_idx] < old_parcel_cloud->r_bubble[p_idx])
    {
        old_parcel_cloud->r_bubble[p_idx] = 0.95 * old_parcel_cloud->radius[p_idx];
    }

    // old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    // printf("running thermal breakup routine p_idx = %i, breakup count = %i \n",p_idx);
    CONVERGE_index_t N = num_child_parcels;

    // End of initialization section
    // printf("\nbreakup count %i",breakup_counter);
    if (old_parcel_cloud->radius[p_idx] < old_parcel_cloud->r_bubble[p_idx])
    {
        old_parcel_cloud->r_bubble[p_idx] = 0.95 * old_parcel_cloud->radius[p_idx];
    }
    // Section 2: Velocity calculations
    
    // Create velocity vectors for all child parcels
    // Get parent parcel's velocity - v = vx i + vy j + vz k
    CONVERGE_precision_t parent_vx, parent_vy, parent_vz;
    CONVERGE_vec3_t parent_velocity,parent_velocity_unit,parent_normal ;
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], &parent_velocity);
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], &parent_velocity_unit); // Parent velocity vector

            // --- Safety check: clamp parent velocity if unphysical ---
        CONVERGE_precision_t parent_vel_mag = sqrt(parent_velocity[0]*parent_velocity[0] +
                                    parent_velocity[1]*parent_velocity[1] +
                                    parent_velocity[2]*parent_velocity[2]);

        // Define a threshold based on your expected flow speeds
        // (e.g. 1e3 m/s is likely way too high for sprays)
        const CONVERGE_precision_t parent_vel_limit = 1.0e3;

        if (parent_vel_mag > parent_vel_limit) {
            printf("Breakup.c WARNING: Parent velocity too large (%e). "
                "Clamping to zero.\n", parent_vel_mag);

            for (int i = 0; i < 3; i++) {
                if(fabs(parent_velocity[i])>parent_vel_limit){
                    parent_velocity[i] = parent_vel_limit;
                }
            }
        }




    CONVERGE_vec3_normalize(parent_velocity_unit);      //Parent Unit velocity vector 
    if(CONVERGE_vec3_length(parent_velocity_unit)<0.99){
        printf("\n Breakup.c\n uu[p_idx] = %e %e %e\nparent_velocity = %e %e %e\nparent_velocity_unit = %e %e %e\n", old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2],parent_velocity[0], parent_velocity[1], parent_velocity[2], parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);
        printf("\n length of parent_velocity_unit = %e",CONVERGE_vec3_length(parent_velocity_unit));
        CONVERGE_mpi_abort();
    }else if(CONVERGE_vec3_length(parent_velocity_unit)>1.01){
        printf("\n Breakup.c\n uu[p_idx] = %e %e %e\nparent_velocity = %e %e %e\nparent_velocity_unit = %e %e %e\n", old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2],parent_velocity[0], parent_velocity[1], parent_velocity[2], parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);
        printf("\n length of parent_velocity_unit = %e",CONVERGE_vec3_length(parent_velocity_unit));
        CONVERGE_mpi_abort();
    }
    // printf("\n Breakup.c\n uu[p_idx] = %e %e %e\nparent_velocity = %e %e %e\nparent_velocity_unit = %e %e %e\n", old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2],parent_velocity[0], parent_velocity[1], parent_velocity[2], parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);


    // CONVERGE_precision_t parent_vmag = CONVERGE_sqrt(CONVERGE_square(parent_vx) + CONVERGE_square(parent_vy) + CONVERGE_square(parent_vz));
    CONVERGE_precision_t parent_vmag = CONVERGE_vec3_length(parent_velocity);
    // Calculate magnitude of child parcel velocity
    CONVERGE_precision_t rad_vel = 3.0 * old_parcel_cloud->v_bubble[p_idx] * CONVERGE_square(old_parcel_cloud->r_bubble[p_idx]) * (old_parcel_cloud->radius[p_idx] - old_parcel_cloud->r_bubble[p_idx]) / (CONVERGE_cube(old_parcel_cloud->radius[p_idx]) - CONVERGE_cube(old_parcel_cloud->r_bubble[p_idx]));
    if (rad_vel > parent_vmag)
    {
        // printf("\nLarge rad vel ---- parent vel magnitude = %e, child_rad_vel = %e", parent_vmag, rad_vel);
        // printf("\n p_idx = %i",p_idx);
        // printf("\n r_bubble = %e",old_parcel_cloud->r_bubble[p_idx]);
        // printf("\n radius = %e",old_parcel_cloud->radius[p_idx]);
        // printf("\n v_bubble = %e",old_parcel_cloud->v_bubble[p_idx]);   
        // printf("\n thermal_breakup_flag = %i",old_parcel_cloud->thermal_breakup_flag[p_idx]);
        // printf("\n tbt = %i",old_parcel_cloud->tbt[p_idx]);
        // printf("\n pbt = %e",old_parcel_cloud->pbt[p_idx]);
        // printf("\n parent velocity = %e %e %e",old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2]);
        // // rad_vel = 0.0;           //Tried this to fix probllem with YNH3 field but it didn't work
    
    }
    else if(fabs(rad_vel)<1.0e-9){
        printf("\n rad_vel = %e",rad_vel);
        printf("\n p_idx = %li",p_idx);
        printf("\n r_bubble = %e",old_parcel_cloud->r_bubble[p_idx]);
        printf("\n radius = %e",old_parcel_cloud->radius[p_idx]);
        printf("\n v_bubble = %e",old_parcel_cloud->v_bubble[p_idx]);   
        printf("\n breakup_phase = %i",old_parcel_cloud->breakup_phase[p_idx]);
        CONVERGE_mpi_abort();
    }


////Start of old child velocity _________________
    // printf("rad _vel =  %e, vmag = %e",rad_vel,parent_vmag);
//perpendicular vector calculation
CONVERGE_vec3_t arbitrary;
do {
    arbitrary[0] = 2.0 * CONVERGE_random_precision() - 1.0;
    arbitrary[1] = 2.0 * CONVERGE_random_precision() - 1.0;
    arbitrary[2] = 2.0 * CONVERGE_random_precision() - 1.0;
    CONVERGE_vec3_normalize(arbitrary);
} while (fabs(CONVERGE_vec3_dot(arbitrary, parent_velocity_unit)) > 0.7);


CONVERGE_vec3_normalize(arbitrary);
CONVERGE_vec3_cross(parent_velocity_unit, arbitrary, &parent_normal);
CONVERGE_vec3_normalize(parent_normal);
CONVERGE_precision_t normal_length = CONVERGE_vec3_length(parent_normal);



// Debug check - should be very close to 1.0
if (fabs(normal_length - 1.0) > 1.0e-1) {
    printf("ERROR: parent_normal not properly normalized! Length = %e\n", normal_length);
    printf("parent_normal = [%e, %e, %e]\n", 
           parent_normal[0], parent_normal[1], parent_normal[2]);
           printf("\n parent_velocity_unit = %e %e %e\n ", parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);
           printf("\n arbritrary = %e %e %e\n ", arbitrary[0], arbitrary[1], arbitrary[2]);
    // CONVERGE_mpi_abort();
}
    // printf("\nparent_normal = %e %e %e\n", parent_normal[0], parent_normal[1], parent_normal[2]);

    
   //-----------------------------Calculate child parcel velocities------------------------------------------------

    // First child parcel will have radial velocity along normal
    CONVERGE_vec3_dup(parent_normal,&user_child_velocity[0]); // Set first child parcel's velocity to be along the normal
    CONVERGE_vec3_normalize(user_child_velocity[0]);
    CONVERGE_vec3_scale(user_child_velocity[0], rad_vel * breakup_velocity_scale);
    


    // Rotate first parcel's velocity by a random angle
    CONVERGE_precision_t rand = CONVERGE_random_precision(); // Generates a random number between 0 and 1
    CONVERGE_precision_t psi = rand * 2 * PI;                // Random angle between 0 and 2 PI
    // printf("psi = %f deg\n",psi*360/(2*PI));
    //  Each of the remaining child parcels will be evenly distributed around the plane perpendicular to the normal
    CONVERGE_precision_t sin_psi = sin(psi);
    CONVERGE_precision_t cos_psi = cos(psi);
    //Rotation around parent's velocity vector by angle psi - Rodrigues' rotation formula
    // child velocity[i] = user_child_velocity[i-1]*cos(psi) + parent_velocity_normal_x_user_child_velocity[i-1] * sin(psi) + parent_velocity_x_parent_velocity_x_user_child_velocity[i-1] * (1 - cos(psi))
    CONVERGE_vec3_t a,b,c,d;

    CONVERGE_vec3_dup(user_child_velocity[0], &a); // Previous child parcel's velocity 

    CONVERGE_vec3_cross(parent_velocity_unit, user_child_velocity[0],&b); 
    CONVERGE_vec3_dup(parent_velocity_unit, &c); // Parent velocity unit vector
    
    CONVERGE_vec3_scale(a, cos_psi) ; //Term 1 
    CONVERGE_vec3_scale(b, sin_psi); // Term 2 
    CONVERGE_vec3_scale(c,CONVERGE_vec3_dot(parent_velocity_unit, user_child_velocity[0])* (1- cos_psi)); //Term 3

    CONVERGE_vec3_add(a,b,&d);
    CONVERGE_vec3_add(d,c, &user_child_velocity[0]); // Final child velocity vector

    CONVERGE_precision_t sin_theta, cos_theta, theta;
    theta = 2 * PI / N; // Angle between child parcels
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    //Developed to split parent parcel into N smaller parcels at breakup 
    for (int jj = 1; jj < N; jj++) // For all the other parcels 2:N
    {
    CONVERGE_vec3_dup(user_child_velocity[jj-1], &a); // Previous child parcel's velocity 

    CONVERGE_vec3_cross(parent_velocity_unit, user_child_velocity[jj-1],&b); 
    CONVERGE_vec3_dup(parent_velocity_unit, &c); // Parent velocity unit vector
    
    CONVERGE_vec3_scale(a, cos_theta) ; //Term 1 
    CONVERGE_vec3_scale(b, sin_theta); // Term 2 
    CONVERGE_vec3_scale(c,CONVERGE_vec3_dot(parent_velocity_unit, user_child_velocity[jj-1])* (1- cos_theta)); //Term 3

    CONVERGE_vec3_add(a,b,&d);
    CONVERGE_vec3_add(d,c, &user_child_velocity[jj]); // Final child velocity vector
    CONVERGE_vec3_normalize(user_child_velocity[jj]);
    if(rad_vel * breakup_velocity_scale > 100.0){
    //  printf("|Vc| = %f",rad_vel * breakup_velocity_scale);   
    }
    CONVERGE_vec3_scale(user_child_velocity[jj],rad_vel * breakup_velocity_scale);
    
    for (int k = 0; k < 3; k++)
    user_child_velocity[jj][k] *= (1.0 + 0.1 * (CONVERGE_random_precision() - 0.5));
   
    } // end of jj loop

  //End of old child initilisation 


  

    //----------------------------Calculate post breakup radius and number of drops for each child parcel----------------------------

    // Radius and Num Drop
    //Pre-Brekaup radius compare
   // printf("\nradius = %e  r_drop_0 = %e r_therm =%e dgre_cycles = %i",old_parcel_cloud->radius[p_idx],old_parcel_cloud->r_drop_0[p_idx],old_parcel_cloud->r_therm[p_idx],old_parcel_cloud->dgre_cycle_count[p_idx]);
    CONVERGE_precision_t rad_denom, rad_term1, rad_term2, rad_term3, rad_term4, parent_radius,parent_nd;
    parent_radius = old_parcel_cloud->radius[p_idx];
    parent_nd= old_parcel_cloud->num_drop[p_idx];
    CONVERGE_precision_t r_bubble_cube = CONVERGE_cube(old_parcel_cloud->r_bubble[p_idx]);
    
    CONVERGE_precision_t denom = CONVERGE_cube(old_parcel_cloud->radius[p_idx]) - r_bubble_cube;
    if (fabs(denom) < 1.0e-20) {
        printf("\nBreakup.c: Error: Denominator in rad_denom calculation is close to zero. Aborting.\n");
        printf("\n radius = %e",old_parcel_cloud->radius[p_idx]);
        printf("\n rububble = %e",old_parcel_cloud->r_bubble[p_idx]);
        CONVERGE_mpi_abort();
    }

    //this is a comment 

    rad_denom = 0.5 * (1.0 / denom);
    rad_term1 = (CONVERGE_square(old_parcel_cloud->radius[p_idx]) + CONVERGE_square(old_parcel_cloud->r_bubble[p_idx]));
    rad_term2 = 3.0 * CONVERGE_square(old_parcel_cloud->v_bubble[p_idx]) * (r_bubble_cube - (r_bubble_cube*old_parcel_cloud->r_bubble[p_idx]) * (1 /parent_radius));
    if (fabs(old_parcel_cloud->surf_ten[p_idx]) < 1.0e-12) {
        printf("\nBreakup.c: Error: Surface tension is close to zero (%e) for parcel %li. Aborting.\n", old_parcel_cloud->surf_ten[p_idx], p_idx);
        CONVERGE_mpi_abort();
    }
    rad_term3 = old_parcel_cloud->density[p_idx] / (3.0 * old_parcel_cloud->surf_ten[p_idx]);

    //printf("\nrad term 3 = %e, den = %e, surten = %e", rad_term3, old_parcel_cloud->density[p_idx], old_parcel_cloud->radius[p_idx]);
    // printf("\nTERM 2 V_BUBBLE = %e R_BUBBLE = %e R_DROP = %e",old_parcel_cloud->v_bubble[p_idx],old_parcel_cloud->r_bubble[p_idx],old_parcel_cloud->radius[p_idx]);
    rad_term4 = CONVERGE_square(breakup_velocity_scale*rad_vel) / 2.0;
    CONVERGE_precision_t B = breakup_radius_scale;           //Constant determining radius of children at breakup (1.0 is default)
CONVERGE_precision_t radius_denominator = 2.0 * B * rad_denom * rad_term1 + rad_term3 * (rad_term2 * rad_denom - rad_term4);
if (fabs(radius_denominator) < 1.0e-20) {
    printf("\nBreakup.c: Error: Denominator for calculated_radius is close to zero (%e) for parcel %li. Aborting.\n", radius_denominator, p_idx);
    CONVERGE_mpi_abort();
}
CONVERGE_precision_t calculated_radius = 1.0 / radius_denominator;
// calculated_radius = calculated_radius * 0.1; // Testing decimating the radius further to see if intensifies evap issues
    // old_parcel_cloud->radius_tm1[p_idx] = old_parcel_cloud->radius[p_idx];
    if (calculated_radius < 0.0)
    {
        printf("\nBreakup.c: Error: calculated_radius is negative. Aborting.\n");
        printf("\nparent radius = %e",parent_radius);
        printf("\n bubble radius = %e",old_parcel_cloud->r_bubble[p_idx]);
        printf("\nrad_denom = %e",rad_denom);
        printf("\nrad_term1 = %e",rad_term1);
        printf("\nrad_term2 = %e",rad_term2);
        printf("\nrad_term3 = %e",rad_term3);
        printf("\nrad_term4 = %e",rad_term4);
        printf("\nrad_vel = %e",rad_vel);
        printf("\nrt3*(rt2*rd -rt4) = %e",rad_term3*(rad_term2*rad_denom - rad_term4));
        CONVERGE_mpi_abort();
    }
    // printf("\ntbt=%i",old_parcel_cloud->tbt[p_idx]);
    if (old_parcel_cloud->r_bubble[p_idx] > parent_radius)
    {
        printf("\nbubble radius larger than droplet's original radius");
    }
    if (calculated_radius > parent_radius *1.01)
    {
        CONVERGE_int_t rankb;
        CONVERGE_mpi_comm_rank(&rankb);
        printf("\n Thermal Breakup has increased radius");
        printf("\n rank = %i p_idx = %li",rankb,p_idx);
        printf("\nr_old = %e r_new = %e rb = %e vb = %e", parent_radius, old_parcel_cloud->radius[p_idx], old_parcel_cloud->r_bubble[p_idx], old_parcel_cloud->v_bubble[p_idx]);
       
        printf("\nrad_term1 = %e rt2 = %e rt3 = %e rt4 = %e rd = %e", rad_term1, rad_term2, rad_term3, rad_term4, rad_denom);
        printf("\nbreakup_phase = %i", old_parcel_cloud->breakup_phase[p_idx]);
        printf("\n vb = %e ",old_parcel_cloud->v_bubble[p_idx]);
        printf("\n Vc = %e rho/sigma = %e", rad_vel, 3.0 * rad_term3);
        printf("\nrt1*rd = %e", rad_term1 * rad_denom);
        printf("\nrt2*rt3*rd = %e", rad_term2 * rad_term3 * rad_denom);
        printf("\nrt4*rt3 = %e", rad_term4 * rad_term3);
        printf("\n DGRE cycle count = %i",old_parcel_cloud->dgre_cycle_count[p_idx]);
        CONVERGE_precision_t P_sat;
        CONVERGE_precision_t Td = old_parcel_cloud->temp[p_idx];
        Saturation_PressureNH3(Td,&P_sat);
        printf("\n P_sat = %e,Td = %f",P_sat,Td);
         printf("\n recalculating v_bubble...");
         CONVERGE_precision_t P_amb = 1.5e6; 
         Bubble_Velocity(old_parcel_cloud,p_idx,P_sat,P_amb);
        printf("\n v_bubble = %e",old_parcel_cloud->v_bubble[p_idx]);
        // printf("\nsurf ten = %e, density = %e rt3 = %e",old_parcel_cloud->surf_ten[p_idx],old_parcel_cloud->density[p_idx],rad_term3);
    }



    prof_calcs += CONVERGE_mpi_wtime() - t0;
    
    // ============================================================================
    // ROSIN-RAMMLER CHILD SAMPLING
    // ============================================================================
    
    // Ensure RR parameters are initialized
    if (!rr_params.initialized) {
        init_RR_distribution(3.2);  // Initialize with default n_RR = 3.2
    }
    
    // Arrays for sampled children
    CONVERGE_precision_t child_radii[MAX_NUM_CHILDREN];
    CONVERGE_precision_t child_num_drop[MAX_NUM_CHILDREN];
    
    // Call RR sampling function with error handling
    int rr_status = sample_RR_children(
        parent_radius,                      // Parent droplet radius
        old_parcel_cloud->num_drop[p_idx],  // Parent num_drop
        calculated_radius,                   // Target R32
        num_child_parcels,                   // N = 12
        rr_params.n_RR,                     // Shape parameter
        rr_params.gamma_ratio,              // Pre-computed gamma ratio
        child_radii,                        // Output: child radii
        child_num_drop                      // Output: child num_drop
    );
    
    // Fallback to uniform distribution if RR sampling failed
    if (rr_status != 0) {
        fallback_uniform_children(
            parent_radius,
            old_parcel_cloud->num_drop[p_idx],
            calculated_radius,
            num_child_parcels,
            child_radii,
            child_num_drop
        );
    }
    
    // ============================================================================
    // END ROSIN-RAMMLER SAMPLING
    // ============================================================================
    
    // Log breakup event data for each child parcel
    static FILE* breakup_events_log2 = NULL;
    if (!breakup_events_log2) {
        breakup_events_log2 = fopen("breakup_events.csv", "a");
    }
    if (breakup_events_log2) {
        for (int i = 0; i < num_child_parcels; i++) {
            fprintf(breakup_events_log2, "%.6e,%ld,%ld,%.6e,%.6e,%.6e,%.6e\n",
                    CONVERGE_simulation_time_sec(), CONVERGE_ncyc(), p_idx,
                    parent_radius, old_parcel_cloud->r_bubble[p_idx], 
                    old_parcel_cloud->v_bubble[p_idx], child_radii[i]);
        }
        fflush(breakup_events_log2);
    }
    
    // DIAGNOSTIC: Calculate mean child radius and log if > 80 μm
    CONVERGE_precision_t mean_child_radius = 0.0;
    for (int i = 0; i < num_child_parcels; i++) {
        mean_child_radius += child_radii[i];
    }
    mean_child_radius /= num_child_parcels;
    
    if (mean_child_radius > 80.0e-6) {
        printf("BREAKUP_LARGE_CHILDREN: t=%.6e, p_idx=%ld, mean_child_R=%.2f um, parent_R=%.2f um\n",
               CONVERGE_simulation_time_sec(), p_idx, mean_child_radius*1e6, old_r*1e6);
        
        static FILE* breakup_log = NULL;
        if (!breakup_log) {
            breakup_log = fopen("breakup_debug.csv", "a");
        }
        if (breakup_log) {
            fprintf(breakup_log, "%.6e,%ld,%ld,LARGE_CHILDREN,%.6e,%.6e,0,0,%.6e,0,%.2f,%d,%.6e\n",
                    CONVERGE_simulation_time_sec(), CONVERGE_ncyc(), p_idx,
                    old_r, old_r_bubble, old_nd, old_parcel_cloud->temp[p_idx],
                    old_parcel_cloud->film_flag[p_idx], mean_child_radius);
            fflush(breakup_log);
        }
    }
    
    //--------- Testing Child Parcel Introduction ----------------//
 
    // Calculate number of child parcels
    CONVERGE_precision_t old_mass, new_mass;
    CONVERGE_index_t nnn;
    CONVERGE_precision_t growth_rate, wave_length, radius_equil;
    CONVERGE_precision_t new_parcel_num_drop, new_parcel_mass, new_radius;
    CONVERGE_vec3_t new_parcel_uu;

    // OLD APPROACH (commented out - now using RR):
    // new_radius = calculated_radius;
    // old_mass = old_parcel_cloud->num_drop[p_idx] * 1.3333 * PI * CONVERGE_cube(old_parcel_cloud->radius[p_idx]);
    // new_mass = old_mass / num_child_parcels;
    // new_parcel_num_drop = new_mass / (1.3333 * PI * CONVERGE_cube(new_radius));
    // new_parcel_num_drop = old_parcel_cloud->num_drop[p_idx];

    //Try cooling the parcel to saturation temp -2 to prevent excessive evap 
    // Just doing this manually for 2 bar for now T -> 252 K 
    // old_parcel_cloud->temp[p_idx] = 252.0;


    t0 = CONVERGE_mpi_wtime();
    growth_rate = 0.0;
    wave_length = 0.0;
    CONVERGE_index_t initial_cloud_size = CONVERGE_cloud_size(cloud);
    // printf("\nInitial cloud size = %i",initial_cloud_size);
    // printf("\nParent parcel radius = %e, num_drop = %e", old_parcel_cloud->radius[p_idx], old_parcel_cloud->num_drop[p_idx]);
    if(initial_cloud_size >0)
    {
        CONVERGE_precision_t nd_before_break = old_parcel_cloud->num_drop[p_idx];

                    //Zero parent drop's radius - this triggers CONVERGE to remove the parent 

                    old_parcel_cloud->radius[p_idx] = 0.0; 
                    old_parcel_cloud->radius_tm1[p_idx] = 0.0;
                    old_parcel_cloud->num_drop[p_idx] = 0.0; 
                    old_parcel_cloud->num_drop_tm1[p_idx] = 0.0; 
                    // old_parcel_cloud->xx[p_idx][0] = 1.0; // This put's the parcel outside of the domain, so it will be removed 
            for(nnn = 0; nnn < num_child_parcels; nnn++)
            {
                // Use sampled RR radius and num_drop for this child
                new_radius = child_radii[nnn];
                new_parcel_num_drop = child_num_drop[nnn];
              
                CONVERGE_vec3_add(parent_velocity, user_child_velocity[nnn], &new_parcel_uu);

                if (CONVERGE_vec3_length(new_parcel_uu) > 1.0e3) {
                    printf("\nBreakup.c: Child velocity too large = %e %e %e. Capping to parent velocity.", new_parcel_uu[0], new_parcel_uu[1], new_parcel_uu[2]);
                    CONVERGE_vec3_dup(parent_velocity, &new_parcel_uu);
                }

                int rank;
	            CONVERGE_mpi_comm_rank(&rank);
      
                user_child_velocity_x = new_parcel_uu[0];
                user_child_velocity_y = new_parcel_uu[1];
                user_child_velocity_z = new_parcel_uu[2];
        

                    CONVERGE_precision_t t1 = CONVERGE_mpi_wtime();
               CONVERGE_spray_child_parcel(new_parcel_uu,
                                            growth_rate,
                                            wave_length,
                                            new_radius,              // Now varies per child (from RR)
                                            new_parcel_num_drop,     // Now from RR sampling
                                            p_idx,
                                            cloud);
                prof_child_parcel += CONVERGE_mpi_wtime() - t1;

            // reload after adding parcels
            load_user_cloud(old_parcel_cloud, cloud);
            }


            CONVERGE_index_t new_cloud_size = CONVERGE_cloud_size(cloud);
            // printf("\nNew cloud size = %i\n\n",new_cloud_size);
            // if(new_cloud_size <= initial_cloud_size)
            // {
            //     printf("\nError: New cloud size is not larger than initial cloud size after breakup");
            //     CONVERGE_mpi_abort();
            // }
        }
            // --------- End of Testing Child Parcel Introduction ----------------//

            prof_loop += CONVERGE_mpi_wtime() - t0;

    // old_parcel_cloud->num_drop[p_idx] = old_parcel_cloud->num_drop[p_idx] * CONVERGE_cube(parent_radius / old_parcel_cloud->radius[p_idx]);
    CONVERGE_precision_t mnew = old_parcel_cloud->num_drop[p_idx] * 1.3333 * PI * CONVERGE_cube(old_parcel_cloud->radius[p_idx]);
    //printf("\nm0 = %e m_old = %e m_new = %e",old_parcel_cloud->m0[p_idx],parent_nd*1.3333*PI*CONVERGE_cube(parent_radius),mnew);
    if(mnew > 0.1*1.01* old_parcel_cloud->m0[p_idx])
    {
        printf("\nBreakup Model has increased droplet mass!!!\n m_old = %e m_new = %e\nnd_old = %e nd_new = %e\nr_old = %e r_new = %e\n Aborting!!!!!!!!\n",old_parcel_cloud->m0[p_idx],mnew,old_nd,old_parcel_cloud->num_drop[p_idx],old_r,old_parcel_cloud->radius[p_idx]);
        CONVERGE_mpi_abort();
    }
    old_parcel_cloud->breakup_phase[p_idx] = 5;  // Mark as child (COMPLETE)
    // old_parcel_cloud->kb[p_idx]=0;
    old_parcel_cloud->int_omega[p_idx]=0.0;
    old_parcel_cloud->r_bubble[p_idx]=0.0;

    if (old_parcel_cloud->num_drop[p_idx] < 0.0)
    {
        printf("\nParcel N_drop < 0!!!!\n");
    }
    // End of child parcel section
    



    int ncyc = CONVERGE_ncyc();
    if (ncyc != last_cycle) {
    int rank;
    CONVERGE_mpi_comm_rank(&rank);
    // printf("Rank %d, Cycle %d Breakup profiling (s): calc=%%f, child_parcel=%%f, loop=%%f\n",
        //    rank, ncyc, prof_loop, prof_child_parcel, prof_calcs);

    double total = prof_loop + prof_calcs;
    // printf("Rank %d, Cycle %d Breakup profiling (%%): calc=%%f, child_parcel=%%f, loop=%%f\n",
        //    rank, ncyc,
        //    100.0*prof_calcs/total,
        //    100.0*prof_child_parcel/total,
        //    100.0*prof_loop/total);

    prof_loop = prof_child_parcel = prof_calcs = 0.0;
    last_cycle = ncyc;
}

    // DIAGNOSTIC: Log completion of breakup for large parcels
    if (old_r > 80.0e-6) {
        CONVERGE_precision_t new_r = old_parcel_cloud->radius[p_idx];
        CONVERGE_precision_t new_r_bubble = old_parcel_cloud->r_bubble[p_idx];
        CONVERGE_precision_t new_nd = old_parcel_cloud->num_drop[p_idx];
        
        // printf("BREAKUP_COMPLETE: p_idx=%ld, R_drop: %.2f->%.2f um, R_bubble: %.2f->%.2f um, num_drop: %.3e->%.3e\n",
        //        p_idx, old_r*1e6, new_r*1e6, old_r_bubble*1e6, new_r_bubble*1e6, old_nd, new_nd);
        
        static FILE* breakup_log = NULL;
        if (!breakup_log) {
            breakup_log = fopen("breakup_debug.csv", "a");
        }
        if (breakup_log) {
            fprintf(breakup_log, "%.6e,%ld,%ld,EXIT_BREAKUP,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.2f,%d\n",
                    CONVERGE_simulation_time_sec(), CONVERGE_ncyc(), p_idx,
                    old_r, old_r_bubble, new_r, new_r_bubble,
                    old_nd, new_nd, old_parcel_cloud->temp[p_idx],
                    old_parcel_cloud->film_flag[p_idx]);
            fflush(breakup_log);
        }
    }
}
