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

/*
 * breakup_phase states:
 *   0 = DISABLED  (parent, not eligible - subcooled, too small, etc.)
 *   1 = ELIGIBLE  (parent, superheated, ready to enter thermal breakup)
 *   2 = ACTIVE    (parent, growing bubble in sub-timestep loop)
 *   3 = RECOVERY  (parent, bubble collapsed, attempting recovery)
 *   4 = READY     (parent, bubble at threshold, ready to fragment)
 *   5 = COMPLETE  (child - result of actual breakup)
 *   6 = BYPASSED  (child - breakup bypassed, reset to injection state)
 */

#include "lagrangian/env.h"

#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <counter.h>
#include <PsatNH3.h>
#include <user_header.h>
#include <globals.h>
#include "Breakup.h"

void print_all_parcel_fields(FILE *fp, const char* parcel_type, int parcel_idx, struct ParcelCloud* parcel_cloud, int num_parcel_species) {
    fprintf(fp, "--- %s Data (Index: %d) ---\n", parcel_type, parcel_idx);
    fprintf(fp, "user_temp_starm1: %e\n", parcel_cloud->user_temp_starm1[parcel_idx]);
    fprintf(fp, "r_bubble: %e\n", parcel_cloud->r_bubble[parcel_idx]);
    fprintf(fp, "v_bubble: %e\n", parcel_cloud->v_bubble[parcel_idx]);
    fprintf(fp, "v_drop: %e\n", parcel_cloud->v_drop[parcel_idx]);
    fprintf(fp, "r_bubble_0: %e\n", parcel_cloud->r_bubble_0[parcel_idx]);
    fprintf(fp, "r_bubble_tm1: %e\n", parcel_cloud->r_bubble_tm1[parcel_idx]);
    fprintf(fp, "v_bubble_tm1: %e\n", parcel_cloud->v_bubble_tm1[parcel_idx]);
    fprintf(fp, "r_drop_0: %e\n", parcel_cloud->r_drop_0[parcel_idx]);
    fprintf(fp, "r_therm: %e\n", parcel_cloud->r_therm[parcel_idx]);
    fprintf(fp, "omega: %e\n", parcel_cloud->omega[parcel_idx]);
    fprintf(fp, "omega_tm1: %e\n", parcel_cloud->omega_tm1[parcel_idx]);
    fprintf(fp, "int_omega: %e\n", parcel_cloud->int_omega[parcel_idx]);
    fprintf(fp, "m0: %e\n", parcel_cloud->m0[parcel_idx]);
    fprintf(fp, "eta_drop: %e\n", parcel_cloud->eta_drop[parcel_idx]);
    fprintf(fp, "eta_drop_0: %e\n", parcel_cloud->eta_drop_0[parcel_idx]);
    fprintf(fp, "user_lag_var_i: %d\n", parcel_cloud->user_lag_var_i[parcel_idx]);
    fprintf(fp, "breakup_phase: %d\n", parcel_cloud->breakup_phase[parcel_idx]);
    fprintf(fp, "child_index: %d\n", parcel_cloud->child_index[parcel_idx]);
    fprintf(fp, "dgre_cycle_count: %d\n", parcel_cloud->dgre_cycle_count[parcel_idx]);
    fprintf(fp, "parcel_index: %d\n", parcel_cloud->parcel_index[parcel_idx]);
    fprintf(fp, "cloud_index: %d\n", parcel_cloud->cloud_index[parcel_idx]);
    fprintf(fp, "user_lag_var_v3: %e, %e, %e\n", parcel_cloud->user_lag_var_v3[parcel_idx][0], parcel_cloud->user_lag_var_v3[parcel_idx][1], parcel_cloud->user_lag_var_v3[parcel_idx][2]);
    fprintf(fp, "user_lag_var_v3b: %e, %e, %e\n", parcel_cloud->user_lag_var_v3b[parcel_idx][0], parcel_cloud->user_lag_var_v3b[parcel_idx][1], parcel_cloud->user_lag_var_v3b[parcel_idx][2]);
    fprintf(fp, "child_uu: %e, %e, %e\n", parcel_cloud->child_uu[parcel_idx][0], parcel_cloud->child_uu[parcel_idx][1], parcel_cloud->child_uu[parcel_idx][2]);
    fprintf(fp, "from_injector: %d\n", parcel_cloud->from_injector[parcel_idx]);
    fprintf(fp, "from_injector_type: %d\n", parcel_cloud->from_injector_type[parcel_idx]);
    fprintf(fp, "on_triangle: %d\n", parcel_cloud->on_triangle[parcel_idx]);
    fprintf(fp, "from_nozzle: %d\n", parcel_cloud->from_nozzle[parcel_idx]);
    fprintf(fp, "film_flag: %d\n", parcel_cloud->film_flag[parcel_idx]);
    fprintf(fp, "just_hit: %d\n", parcel_cloud->just_hit[parcel_idx]);
    fprintf(fp, "just_hit_leiden: %d\n", parcel_cloud->just_hit_leiden[parcel_idx]);
    fprintf(fp, "is_thick: %c\n", parcel_cloud->is_thick[parcel_idx]);
    fprintf(fp, "rel_vel: %e, %e, %e\n", parcel_cloud->rel_vel[parcel_idx][0], parcel_cloud->rel_vel[parcel_idx][1], parcel_cloud->rel_vel[parcel_idx][2]);
    fprintf(fp, "uprime: %e, %e, %e\n", parcel_cloud->uprime[parcel_idx][0], parcel_cloud->uprime[parcel_idx][1], parcel_cloud->uprime[parcel_idx][2]);
    fprintf(fp, "uu: %e, %e, %e\n", parcel_cloud->uu[parcel_idx][0], parcel_cloud->uu[parcel_idx][1], parcel_cloud->uu[parcel_idx][2]);
    fprintf(fp, "uu_tm1: %e, %e, %e\n", parcel_cloud->uu_tm1[parcel_idx][0], parcel_cloud->uu_tm1[parcel_idx][1], parcel_cloud->uu_tm1[parcel_idx][2]);
    fprintf(fp, "xx: %e, %e, %e\n", parcel_cloud->xx[parcel_idx][0], parcel_cloud->xx[parcel_idx][1], parcel_cloud->xx[parcel_idx][2]);
    fprintf(fp, "xx_tm1: %e, %e, %e\n", parcel_cloud->xx_tm1[parcel_idx][0], parcel_cloud->xx_tm1[parcel_idx][1], parcel_cloud->xx_tm1[parcel_idx][2]);
    fprintf(fp, "v_nu: %e\n", parcel_cloud->v_nu[parcel_idx]);
    fprintf(fp, "v_sh: %e\n", parcel_cloud->v_sh[parcel_idx]);
    fprintf(fp, "temp: %e\n", parcel_cloud->temp[parcel_idx]);
    fprintf(fp, "temp_tm1: %e\n", parcel_cloud->temp_tm1[parcel_idx]);
    fprintf(fp, "temp_starm1: %e\n", parcel_cloud->temp_starm1[parcel_idx]);
    fprintf(fp, "rey_num: %e\n", parcel_cloud->rey_num[parcel_idx]);
    fprintf(fp, "rel_vel_mag: %e\n", parcel_cloud->rel_vel_mag[parcel_idx]);
    fprintf(fp, "radius: %e\n", parcel_cloud->radius[parcel_idx]);
    fprintf(fp, "radius_tm1: %e\n", parcel_cloud->radius_tm1[parcel_idx]);
    fprintf(fp, "parent: %e\n", parcel_cloud->parent[parcel_idx]);
    fprintf(fp, "density: %e\n", parcel_cloud->density[parcel_idx]);
    fprintf(fp, "density_tm1: %e\n", parcel_cloud->density_tm1[parcel_idx]);
    fprintf(fp, "gas_density: %e\n", parcel_cloud->gas_density[parcel_idx]);
    for (int i = 0; i < num_parcel_species; ++i) {
        fprintf(fp, "mfrac[%d]: %e\n", i, parcel_cloud->mfrac[parcel_idx * num_parcel_species + i]);
        fprintf(fp, "mfrac_tm1[%d]: %e\n", i, parcel_cloud->mfrac_tm1[parcel_idx * num_parcel_species + i]);
        fprintf(fp, "dm_dt[%d]: %e\n", i, parcel_cloud->dm_dt[parcel_idx * num_parcel_species + i]);
        fprintf(fp, "drdt[%d]: %e\n", i, parcel_cloud->drdt[parcel_idx * num_parcel_species + i]);
    }
    fprintf(fp, "num_drop: %e\n", parcel_cloud->num_drop[parcel_idx]);
    fprintf(fp, "surf_temp: %e\n", parcel_cloud->surf_temp[parcel_idx]);
    fprintf(fp, "tbreak_kh: %e\n", parcel_cloud->tbreak_kh[parcel_idx]);
    fprintf(fp, "shed_num_drop: %e\n", parcel_cloud->shed_num_drop[parcel_idx]);
    fprintf(fp, "shed_mass: %e\n", parcel_cloud->shed_mass[parcel_idx]);
    fprintf(fp, "sactive: %e\n", parcel_cloud->sactive[parcel_idx]);
    fprintf(fp, "sactive_tm1: %e\n", parcel_cloud->sactive_tm1[parcel_idx]);
    fprintf(fp, "surf_ten: %e\n", parcel_cloud->surf_ten[parcel_idx]);
    fprintf(fp, "viscosity: %e\n", parcel_cloud->viscosity[parcel_idx]);
    fprintf(fp, "distort: %e\n", parcel_cloud->distort[parcel_idx]);
    fprintf(fp, "distort_dot: %e\n", parcel_cloud->distort_dot[parcel_idx]);
    fprintf(fp, "wall_heat_exchange: %e\n", parcel_cloud->wall_heat_exchange[parcel_idx]);
    fprintf(fp, "l_rr: %e\n", parcel_cloud->l_rr[parcel_idx]);
    fprintf(fp, "l_rc: %e\n", parcel_cloud->l_rc[parcel_idx]);
    fprintf(fp, "l_temp1: %e\n", parcel_cloud->l_temp1[parcel_idx]);
    fprintf(fp, "l_temp2: %e\n", parcel_cloud->l_temp2[parcel_idx]);
    fprintf(fp, "film_shed: %e\n", parcel_cloud->film_shed[parcel_idx]);
    fprintf(fp, "tbreak_rt: %e\n", parcel_cloud->tbreak_rt[parcel_idx]);
    fprintf(fp, "surf_temp_tm1: %e\n", parcel_cloud->surf_temp_tm1[parcel_idx]);
    fprintf(fp, "num_drop_tm1: %e\n", parcel_cloud->num_drop_tm1[parcel_idx]);
    fprintf(fp, "film_energy: %e\n", parcel_cloud->film_energy[parcel_idx]);
    fprintf(fp, "t_turb: %e\n", parcel_cloud->t_turb[parcel_idx]);
    fprintf(fp, "t_turb_accum: %e\n", parcel_cloud->t_turb_accum[parcel_idx]);
    fprintf(fp, "film_thickness: %e\n", parcel_cloud->film_thickness[parcel_idx]);
    fprintf(fp, "area_reduction: %e\n", parcel_cloud->area_reduction[parcel_idx]);
    fprintf(fp, "tke0: %e\n", parcel_cloud->tke0[parcel_idx]);
    fprintf(fp, "eps0: %e\n", parcel_cloud->eps0[parcel_idx]);
    fprintf(fp, "lifetime: %e\n", parcel_cloud->lifetime[parcel_idx]);
    if (parcel_cloud->film_accum_bit_flag) fprintf(fp, "film_accum_bit_flag: %llu\n", parcel_cloud->film_accum_bit_flag[parcel_idx]);
    if (parcel_cloud->film_accum_plus_bit_flag) fprintf(fp, "film_accum_plus_bit_flag: %u\n", parcel_cloud->film_accum_plus_bit_flag[parcel_idx]);
    if (parcel_cloud->film_thickness_tm1) fprintf(fp, "film_thickness_tm1: %e\n", parcel_cloud->film_thickness_tm1[parcel_idx]);
    if (parcel_cloud->area_in_film) fprintf(fp, "area_in_film: %e\n", parcel_cloud->area_in_film[parcel_idx]);
}

/**********************************************************************/ 
/*                                                                    */ 
/* Name: user_parcel_prop                                             */ 
/*                                                                    */ 
/* Description: user_parcel_prop sets values for the     custom       */ 
/* parcel properties when new parcels are created. This is            */ 
/* accomplished through the four functions: user_parcel_inject,       */ 
/* user_parcel_child, user_parcel_film, user_parcel_splash            */ 
/*                                                                    */ 
/*                                                                    */ 
/**********************************************************************/ 

// initialize values for the custom parcel properties when new parcels 
// are first injected into the domain
CONVERGE_UDF(parcel_inject,
             IN(VALUE(CONVERGE_index_t, passed_parcel_idx), VALUE(CONVERGE_cloud_t, passed_spray_cloud)),
             OUT(CONVERGE_VOID))
{
   
   // NOTE: Cannot access CFD mesh pressure here (parcel_inject called before drop_distort)
   // Use a reasonable default - will be corrected in RPE solver on first call
   CONVERGE_precision_t ambient_pres = 2.0e5;  // Default 2 bar, recalculated in RPE solver

   // printf("\n START OF PARCEL_PROP.C \n");
   struct ParcelCloud parcel_cloud;

   // printf("PARCEL_PROP.C L58\n");
   load_user_cloud(&parcel_cloud, passed_spray_cloud);
   // Calculate Saturation Pressure from Antoine's Equation
   CONVERGE_precision_t P_sat,Td;
   Td = parcel_cloud.temp[passed_parcel_idx];
   Saturation_PressureNH3(Td, &P_sat);

   // printf("PARCEL_PROP.C L70 P_sat = %f\n",P_sat);
   // printf("PARCEL_PROP.C L67\n");

   // printf("\n ambient_pres = %f \n",ambient_pres);
   // printf("\n r_bubble old = %f \n", parcel_cloud.r_bubble[passed_parcel_idx]);

   // Initialize bubble to 0.0 - RPE solvers will calculate proper R0 with actual P_amb
   // This avoids using default ambient_pres which may not match local CFD pressure
   parcel_cloud.r_bubble[passed_parcel_idx] = 0.0;
   parcel_cloud.r_bubble_0[passed_parcel_idx] = 0.0;

   parcel_cloud.time_of_injection[passed_parcel_idx] = CONVERGE_simulation_time_sec();
   // printf("\n PARCEL_PROP.C L69 r_bubble = %e\n", parcel_cloud.r_bubble[passed_parcel_idx]);
   parcel_cloud.v_bubble[passed_parcel_idx] = 0.0;

   // printf("\n\n r_bubble = %e <seg_60> r_bubble_0 = %e", parcel_cloud.r_bubble[passed_parcel_idx], parcel_cloud.r_bubble_0[passed_parcel_idx]);

   // DIAGNOSTIC: Hijack film_thickness for r_bubble output (safe if not using film model)
   parcel_cloud.film_thickness[passed_parcel_idx] = parcel_cloud.r_bubble[passed_parcel_idx];

   // R_D_0
   parcel_cloud.r_drop_0[passed_parcel_idx] = parcel_cloud.radius[passed_parcel_idx];
   parcel_cloud.temp_drop_0[passed_parcel_idx] = parcel_cloud.temp[passed_parcel_idx];
   // Droplet radius for thermal breakup purposes
   parcel_cloud.r_therm[passed_parcel_idx] = parcel_cloud.radius[passed_parcel_idx];
   // Zero Omega and Eta (breakup variables)
   parcel_cloud.int_omega[passed_parcel_idx] = 0;
   parcel_cloud.omega[passed_parcel_idx] = parcel_cloud.omega_tm1[passed_parcel_idx] = 0;
   parcel_cloud.eta_drop[passed_parcel_idx] = 0;

   //Initialze eta_0 according to log-normal distribution
      // mean target for eta_0
      CONVERGE_precision_t m = 0.05; 

      // choose spread
      CONVERGE_precision_t sigma = 0.1; // adjust 0.1–0.5 as needed

      // convert to log-space parameters
      CONVERGE_precision_t mu = log(m) - 0.5 * sigma * sigma;

      // generate two uniforms
      CONVERGE_precision_t u1 = CONVERGE_random_precision();
      CONVERGE_precision_t u2 = CONVERGE_random_precision();

      // transform to normal(0,1)
      CONVERGE_precision_t z = sqrt(-2.0*log(u1)) * cos(2.0*M_PI*u2);

      // scale/shift to normal(mu, sigma^2)
      CONVERGE_precision_t y = mu + sigma * z;

      // exponentiate to get lognormal
      parcel_cloud.eta_drop_0[passed_parcel_idx] = exp(y) * parcel_cloud.radius[passed_parcel_idx];
      // printf("\n eta_drop_0 = %e\n",parcel_cloud.eta_drop_0[passed_parcel_idx]);


   parcel_cloud.dgre_cycle_count[passed_parcel_idx] = 0;
   parcel_cloud.breakup_phase[passed_parcel_idx] = 1;  // ELIGIBLE - ready to enter thermal breakup
   parcel_cloud.film_flag[passed_parcel_idx] = parcel_cloud.breakup_phase[passed_parcel_idx];  // Hijack: mirror breakup_phase
   parcel_cloud.m0[passed_parcel_idx] = (1.33333 * PI * CONVERGE_cube(parcel_cloud.radius[passed_parcel_idx]) * parcel_cloud.num_drop[passed_parcel_idx]);
   parcel_cloud.parcel_index[passed_parcel_idx] = user_parcel_counter;
   parcel_cloud.cloud_index[passed_parcel_idx] = -1;
   
   // Initialize collapse counter to 0
   parcel_cloud.user_lag_var_i[passed_parcel_idx] = 0;
   
   // Initialize recovery fields
   parcel_cloud.recovery_time[passed_parcel_idx] = 0.0;
   parcel_cloud.recovery_count[passed_parcel_idx] = 0;
   
   user_parcel_counter ++;
   

   char *filename = "DGRE.txt";
   char *filename1 = "Temp_Tracker.txt";
   CONVERGE_index_t ncyc = CONVERGE_ncyc();
   // Populate Text File With Header
   if (ncyc == 12)
   {
      FILE *fp = fopen("DGRE.txt", "w");
      if (fp == NULL)
      {
         printf("Error opening the file %s", filename);
      }
      fprintf(fp, "a   b       c       d       real(x0)      real(x1)      real(x2)      omega   eta\tT_drop\tP_sat\tRb\tRd     kb\tt_parcel\tbreakup_flag\n");
      fclose(fp);
      FILE *fp1 = fopen("Temp_Tracker.txt", "w");
      if (fp1 == NULL)
      {
         printf("Error opening the file %s", filename1);
      }
      fprintf(fp1, "CI    PI    Td    Rd    age    Vmag\n");
      fclose(fp1);
   }
}

// initialize values for the custom parcel properties when new child parcels
// are created from the Kelvin-Helmholtz stripping breakup mechanism

// NOTE: use the following syntax to set the child parcel value equal to its parent's value

CONVERGE_UDF(parcel_child,
             IN(VALUE(CONVERGE_index_t, passed_child_parcel_idx),
                VALUE(CONVERGE_index_t, passed_parent_parcel_idx),
               //  VALUE(CONVERGE_vec3_t, user_child_velocity),
                VALUE(CONVERGE_cloud_t, passed_spray_cloud)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud parcel_cloud;
   CONVERGE_precision_t parcel_semi_mass_old, parcel_semi_mass_new;

   load_user_cloud(&parcel_cloud, passed_spray_cloud);

      //Bookkeeping for Temp_Tracker.txt
      parcel_cloud.parcel_index[passed_child_parcel_idx] = passed_child_parcel_idx;
       parcel_cloud.cloud_index[passed_child_parcel_idx] = parcel_cloud.cloud_index[passed_parent_parcel_idx];

   //In case of KH-RT breakup
   // R_D_0 - should not be used for child parcels
   parcel_cloud.r_drop_0[passed_child_parcel_idx] = parcel_cloud.radius[passed_parent_parcel_idx];
   parcel_cloud.r_therm[passed_child_parcel_idx] = parcel_cloud.radius[passed_parent_parcel_idx];



   
   if (parcel_cloud.breakup_phase[passed_parent_parcel_idx] != 0)  // Parent underwent some form of breakup
   {

      //Zero Breakup Properties
      parcel_cloud.r_bubble[passed_child_parcel_idx] = 0.0;
      parcel_cloud.r_bubble_0[passed_child_parcel_idx] = 0.0;
      parcel_cloud.v_bubble[passed_child_parcel_idx] = 0.0;
      parcel_cloud.r_bubble_0[passed_child_parcel_idx] = 0.0;
      
      // Mark as child (COMPLETE) - will never enter thermal breakup
      parcel_cloud.breakup_phase[passed_child_parcel_idx] = 5;
      parcel_cloud.film_flag[passed_child_parcel_idx] = 5;  // Hijack: mirror breakup_phase
      
      parcel_cloud.temp[passed_child_parcel_idx] = parcel_cloud.temp[passed_parent_parcel_idx];
      parcel_cloud.temp_tm1[passed_child_parcel_idx] = parcel_cloud.temp_tm1[passed_parent_parcel_idx];

      // DIAGNOSTIC: Hijack film_thickness for r_bubble output (safe if not using film model)
      parcel_cloud.film_thickness[passed_child_parcel_idx] = parcel_cloud.r_bubble[passed_child_parcel_idx];

      // R_D_0 - should not be used for child parcels
      parcel_cloud.r_drop_0[passed_child_parcel_idx] = parcel_cloud.radius[passed_child_parcel_idx];
      // Droplet radius for thermal breakup purposes
      parcel_cloud.r_therm[passed_child_parcel_idx] = parcel_cloud.radius[passed_child_parcel_idx];
      // Zero Omega and Eta (breakup variables)
      parcel_cloud.int_omega[passed_child_parcel_idx] = 0;
      parcel_cloud.omega[passed_child_parcel_idx] = parcel_cloud.omega_tm1[passed_child_parcel_idx] = 0;
      parcel_cloud.eta_drop[passed_child_parcel_idx] = 0;
      parcel_cloud.dgre_cycle_count[passed_child_parcel_idx] = 0;
      parcel_cloud.radius_tm1[passed_child_parcel_idx] = parcel_cloud.radius[passed_child_parcel_idx];
      // In the parcel_child function, add this line where other parcel properties are being set:
      parcel_cloud.lifetime[passed_child_parcel_idx] = 0.0;  // Reset lifetime for new child droplets 
      
      // Initialize collapse counter to 0
      parcel_cloud.user_lag_var_i[passed_child_parcel_idx] = 0;
      
      // Initialize recovery fields
      parcel_cloud.recovery_time[passed_child_parcel_idx] = 0.0;
      parcel_cloud.recovery_count[passed_child_parcel_idx] = 0;
      
      // If parent's thermal_breakup_flag is set, displace the child parcel
      int rank;
      CONVERGE_mpi_comm_rank(&rank);

      parcel_cloud.uu[passed_child_parcel_idx][0] = user_child_velocity_x;
      parcel_cloud.uu[passed_child_parcel_idx][1] = user_child_velocity_y;
      parcel_cloud.uu[passed_child_parcel_idx][2] = user_child_velocity_z;
      
      parcel_cloud.uu_tm1[passed_child_parcel_idx][0] = user_child_velocity_x;
      parcel_cloud.uu_tm1[passed_child_parcel_idx][1] = user_child_velocity_y;
      parcel_cloud.uu_tm1[passed_child_parcel_idx][2] = user_child_velocity_z;

    
      // Update positions based on normalized velocity direction
      // Displace child by parent_radius in the direction of child velocity
      CONVERGE_precision_t vel_mag = sqrt(user_child_velocity_x * user_child_velocity_x +
                                           user_child_velocity_y * user_child_velocity_y +
                                           user_child_velocity_z * user_child_velocity_z);
      
      // Only apply displacement if velocity is non-zero
      if (vel_mag > 1.0e-10) {
         // Normalize velocity to get unit vector
         CONVERGE_precision_t unit_x = user_child_velocity_x / vel_mag;
         CONVERGE_precision_t unit_y = user_child_velocity_y / vel_mag;
         CONVERGE_precision_t unit_z = user_child_velocity_z / vel_mag;
         
         // Displace by parent radius in the direction of normalized velocity
         CONVERGE_precision_t displacement = parcel_cloud.radius[passed_parent_parcel_idx];
         parcel_cloud.xx[passed_child_parcel_idx][0] = parcel_cloud.xx[passed_parent_parcel_idx][0] + displacement * unit_x;
         parcel_cloud.xx[passed_child_parcel_idx][1] = parcel_cloud.xx[passed_parent_parcel_idx][1] + displacement * unit_y;
         parcel_cloud.xx[passed_child_parcel_idx][2] = parcel_cloud.xx[passed_parent_parcel_idx][2] + displacement * unit_z;
      } else {
         // Zero velocity - place at same location as parent
         parcel_cloud.xx[passed_child_parcel_idx][0] = parcel_cloud.xx[passed_parent_parcel_idx][0];
         parcel_cloud.xx[passed_child_parcel_idx][1] = parcel_cloud.xx[passed_parent_parcel_idx][1];
         parcel_cloud.xx[passed_child_parcel_idx][2] = parcel_cloud.xx[passed_parent_parcel_idx][2];
      }
      parcel_cloud.xx_tm1[passed_child_parcel_idx][0] = parcel_cloud.xx[passed_child_parcel_idx][0];
      parcel_cloud.xx_tm1[passed_child_parcel_idx][1] = parcel_cloud.xx[passed_child_parcel_idx][1];
      parcel_cloud.xx_tm1[passed_child_parcel_idx][2] = parcel_cloud.xx[passed_child_parcel_idx][2];

   const int num_parcel_species = 1;

   static int child_print_counter = 0;
   if (child_print_counter % 1000 == 0) {
       FILE *fp = fopen("parcel_child_debug.txt", "a");
       if (fp == NULL) {
           printf("Error opening parcel_child_debug.txt\n");
       } else {
           fprintf(fp, "\n\nNCYC: %ld, Time: %e\n", (long)CONVERGE_ncyc(), CONVERGE_simulation_time());
           print_all_parcel_fields(fp, "Parent", passed_parent_parcel_idx, &parcel_cloud, num_parcel_species);
           print_all_parcel_fields(fp, "Child", passed_child_parcel_idx, &parcel_cloud, num_parcel_species);
           fclose(fp);
       }
   }
   child_print_counter++;
   }
   // }
      // parcel_cloud.tbt[passed_child_parcel_idx] = 0;


   parcel_cloud.m0[passed_child_parcel_idx] = (1.33333 * PI * CONVERGE_cube(parcel_cloud.radius[passed_parent_parcel_idx]) * parcel_cloud.num_drop[passed_parent_parcel_idx]);
   //Don't reset breakup flag
   // // printf("\n\n parcel_child called for parcel %i with parent %i which has temp %f\nradius %e and num_drop %e\nparent_radius = %e, parent_num_drop = %e\n",
   //        passed_child_parcel_idx, passed_parent_parcel_idx, parcel_cloud.temp[passed_child_parcel_idx],
   //        parcel_cloud.radius[passed_child_parcel_idx], parcel_cloud.num_drop[passed_child_parcel_idx],
   //        parcel_cloud.radius[passed_parent_parcel_idx], parcel_cloud.num_drop[passed_parent_parcel_idx]);

   // --- GEMINI EDIT: Add validation check for newly created child parcel ---
   int child_phase = parcel_cloud.breakup_phase[passed_child_parcel_idx];
   // A new child should have phase=5 (COMPLETE). 
   // Any other value indicates memory corruption or logic error.
   if (child_phase != 5) {
       CONVERGE_logger_fatal("MEMORY CORRUPTION DETECTED in parcel_prop.c! New child (idx: %d) created with invalid breakup_phase: %d (expected 5)",
                              passed_child_parcel_idx, child_phase);
       
       // Also print the parent's data, as it might provide clues.
       CONVERGE_logger_fatal("Parent (idx: %d) has breakup_phase: %d, radius: %e, temp: %f",
                              passed_parent_parcel_idx, parcel_cloud.breakup_phase[passed_parent_parcel_idx],
                              parcel_cloud.radius[passed_parent_parcel_idx], parcel_cloud.temp[passed_parent_parcel_idx]);

       CONVERGE_mpi_abort();
   }
   // --- END GEMINI EDIT ---
}
// set values for the custom parcel properties when film parcels separate
// from a surface and are converted to spray parcels

// NOTE: set the properties equal to themselves if they are to retain their
//       pre-separation values

CONVERGE_UDF(parcel_splash,
             IN(VALUE(CONVERGE_index_t, passed_spray_parcel_idx),
                VALUE(CONVERGE_index_t, passed_film_parcel_idx),
                VALUE(CONVERGE_cloud_t, passed_spray_cloud),
                VALUE(CONVERGE_cloud_t, passed_film_cloud)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud spray_parcel_cloud;
   struct ParcelCloud film_parcel_cloud;
   CONVERGE_precision_t parcel_semi_mass_old, parcel_semi_mass_new;

   load_user_cloud(&spray_parcel_cloud, passed_spray_cloud);
   load_user_cloud(&film_parcel_cloud, passed_film_cloud);

   parcel_semi_mass_old =
       spray_parcel_cloud.density[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx] *
       spray_parcel_cloud.radius[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx];

   spray_parcel_cloud.density[passed_spray_parcel_idx] = film_parcel_cloud.density[passed_film_parcel_idx];
   spray_parcel_cloud.radius[passed_spray_parcel_idx] = film_parcel_cloud.radius[passed_film_parcel_idx];
   spray_parcel_cloud.density_tm1[passed_spray_parcel_idx] = film_parcel_cloud.density_tm1[passed_film_parcel_idx];

   parcel_semi_mass_new =
       spray_parcel_cloud.density[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx] *
       spray_parcel_cloud.radius[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx];

   spray_parcel_cloud.num_drop[passed_spray_parcel_idx] =
   spray_parcel_cloud.num_drop[passed_spray_parcel_idx] * parcel_semi_mass_old / parcel_semi_mass_new;

      spray_parcel_cloud.parcel_index[passed_spray_parcel_idx] = spray_parcel_cloud.parcel_index[passed_film_parcel_idx];
      spray_parcel_cloud.cloud_index[passed_spray_parcel_idx] = spray_parcel_cloud.cloud_index[passed_film_parcel_idx];
}
// initialize values for the custom parcel properties when new parcels are created
// from film stripping

// NOTE: use the following syntax to set the splashed parcel value equal to the impinged parcel's value

CONVERGE_UDF(parcel_strip,
             IN(VALUE(CONVERGE_index_t, passed_spray_parcel_idx),
                VALUE(CONVERGE_index_t, passed_film_parcel_idx),
                VALUE(CONVERGE_cloud_t, passed_spray_cloud),
                VALUE(CONVERGE_cloud_t, passed_film_cloud)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud spray_parcel_cloud;
   struct ParcelCloud film_parcel_cloud;
   CONVERGE_precision_t parcel_semi_mass_old, parcel_semi_mass_new;

   load_user_cloud(&spray_parcel_cloud, passed_spray_cloud);
   load_user_cloud(&film_parcel_cloud, passed_film_cloud);

   parcel_semi_mass_old =
       spray_parcel_cloud.density[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx] *
       spray_parcel_cloud.radius[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx];

   spray_parcel_cloud.density[passed_spray_parcel_idx] = film_parcel_cloud.density[passed_film_parcel_idx];
   spray_parcel_cloud.radius[passed_spray_parcel_idx] = film_parcel_cloud.radius[passed_film_parcel_idx];
   spray_parcel_cloud.density_tm1[passed_spray_parcel_idx] = film_parcel_cloud.density_tm1[passed_film_parcel_idx];

   parcel_semi_mass_new =
       spray_parcel_cloud.density[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx] *
       spray_parcel_cloud.radius[passed_spray_parcel_idx] * spray_parcel_cloud.radius[passed_spray_parcel_idx];

   spray_parcel_cloud.num_drop[passed_spray_parcel_idx] =
       spray_parcel_cloud.num_drop[passed_spray_parcel_idx] * parcel_semi_mass_old / parcel_semi_mass_new;

         spray_parcel_cloud.parcel_index[passed_spray_parcel_idx] = spray_parcel_cloud.parcel_index[passed_film_parcel_idx];
       spray_parcel_cloud.cloud_index[passed_spray_parcel_idx] = spray_parcel_cloud.cloud_index[passed_film_parcel_idx];
}
