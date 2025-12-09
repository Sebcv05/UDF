/*******************************************************************************
 * CONVERGENT SCIENCE CONFIDENTIAL                                              *
 * All rights reserved.                                                         *
 * All information contained herein is the property of Convergent Science.      *
 * The intellectual and technical concepts contained herein are                 *
 * proprietary to Convergent Science.                                           *
 * Dissemination of this information or reproduction of this material           *
 * is strictly forbidden unless prior written permission is obtained from       *
 * Convergent Science.                                                          *
 *******************************************************************************/

#include <parcel_reset.h>
#include <math.h>
#include <stdio.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

void reset_parcel_to_child(struct ParcelCloud* parcel_cloud, 
                            CONVERGE_index_t p_idx, 
                            const char* reason)
{
    // Get current state
    CONVERGE_precision_t R_old = parcel_cloud->radius[p_idx];
    CONVERGE_precision_t R_new = parcel_cloud->r_drop_0[p_idx];  // Injection radius
    CONVERGE_precision_t Nd_old = parcel_cloud->num_drop[p_idx];
    
    // Calculate total mass (must be conserved)
    // Mass = Nd * (4/3) * PI * R^3 * rho
    // Since rho is constant, we need: Nd_old * R_old^3 = Nd_new * R_new^3
    CONVERGE_precision_t volume_ratio = (R_old * R_old * R_old) / (R_new * R_new * R_new);
    CONVERGE_precision_t Nd_new = Nd_old * volume_ratio;
    
    // Safety check: if injection radius is invalid, use current radius
    if (R_new <= 0.0 || R_new > R_old * 2.0) {
        static int warning_count = 0;
        if (warning_count < 5) {
            printf("[PARCEL_RESET_WARNING] Invalid r_drop_0=%.3e m, using current radius %.3e m (p_idx=%li, reason: %s)\n",
                   R_new, R_old, p_idx, reason);
            warning_count++;
        }
        R_new = R_old;
        Nd_new = Nd_old;
    }
    
    // Perform reset
    parcel_cloud->radius[p_idx] = R_new;
    parcel_cloud->num_drop[p_idx] = Nd_new;
    parcel_cloud->r_bubble[p_idx] = 0.0;
    parcel_cloud->r_bubble_0[p_idx] = 0.0;
    parcel_cloud->v_bubble[p_idx] = 0.0;
    parcel_cloud->breakup_phase[p_idx] = 5;  // Mark as child (COMPLETE)
    parcel_cloud->film_flag[p_idx] = 5;  // Hijack: mirror breakup_phase
    
    // Optional: log the reset (can be disabled if too verbose)
    static int reset_count = 0;
    if (reset_count < 10) {
        printf("[PARCEL_RESET] p_idx=%li, reason: %s\n", p_idx, reason);
        printf("               R_drop: %.3e -> %.3e m (injection radius)\n", R_old, R_new);
        printf("               num_drop: %.3e -> %.3e (mass conserved, ratio=%.3f)\n", 
               Nd_old, Nd_new, volume_ratio);
        printf("               r_bubble: %.3e -> 0.00e+00 m (zeroed)\n", 
               parcel_cloud->r_bubble[p_idx]);
        printf("               breakup_phase: -> 5 (child/complete)\n");
        reset_count++;
    }
}
