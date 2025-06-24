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
#include "lagrangian/env.h"

#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <counter.h>
#include <PsatNH3.h>
#include <user_header.h>

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

   parcel_cloud.r_bubble[passed_parcel_idx] = 2.0 * parcel_cloud.surf_ten[passed_parcel_idx] / (P_sat - ambient_pres);
   if (parcel_cloud.r_bubble[passed_parcel_idx] < 0.0)
   {
      parcel_cloud.r_bubble[passed_parcel_idx] = 0.0;
      parcel_cloud.r_bubble_0[passed_parcel_idx] = 0.0;
   }

   // printf("\n PARCEL_PROP.C L69 r_bubble = %e\n", parcel_cloud.r_bubble[passed_parcel_idx]);
   parcel_cloud.v_bubble[passed_parcel_idx] = 0.0;
   parcel_cloud.r_bubble_0[passed_parcel_idx] = parcel_cloud.r_bubble[passed_parcel_idx];
   // printf("\n END OF PARCEL_PROP.C \n");
   // printf("\n\n r_bubble = %e 	r_bubble_0 = %e", parcel_cloud.r_bubble[passed_parcel_idx], parcel_cloud.r_bubble_0[passed_parcel_idx]);

   // R_D_0
   parcel_cloud.r_drop_0[passed_parcel_idx] = parcel_cloud.radius[passed_parcel_idx];

   // Droplet radius for thermal breakup purposes
   parcel_cloud.r_therm[passed_parcel_idx] = parcel_cloud.radius[passed_parcel_idx];
   // Zero Omega and Eta (breakup variables)
   parcel_cloud.int_omega[passed_parcel_idx] = 0;
   parcel_cloud.omega[passed_parcel_idx] = parcel_cloud.omega_tm1[passed_parcel_idx] = 0;
   parcel_cloud.eta_drop[passed_parcel_idx] = 0;
   parcel_cloud.dgre_cycle_count[passed_parcel_idx] = 0;
   parcel_cloud.tbt[passed_parcel_idx] = 0;
   parcel_cloud.pbt[passed_parcel_idx] = 1;
   parcel_cloud.m0[passed_parcel_idx] = (1.33333 * PI * CONVERGE_cube(parcel_cloud.radius[passed_parcel_idx]) * parcel_cloud.num_drop[passed_parcel_idx]);
   // Set breakup flag to 0
   parcel_cloud.thermal_breakup_flag[passed_parcel_idx] = -1;
   parcel_cloud.parcel_index[passed_parcel_idx] = user_parcel_counter;
   parcel_cloud.cloud_index[passed_parcel_idx] = -1;
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
      fprintf(fp, "a   b       c       d       real(x0)      real(x1)      real(x2)      omega   eta	T_drop	P_sat	Rb	Rd     kb	t_parcel	breakup_flag\n");
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
                VALUE(CONVERGE_cloud_t, passed_spray_cloud)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud parcel_cloud;
   CONVERGE_precision_t parcel_semi_mass_old, parcel_semi_mass_new;

   load_user_cloud(&parcel_cloud, passed_spray_cloud);

   // parcel_semi_mass_old = parcel_cloud.density[passed_child_parcel_idx] * parcel_cloud.radius[passed_child_parcel_idx] *
   //                        parcel_cloud.radius[passed_child_parcel_idx] * parcel_cloud.radius[passed_child_parcel_idx];

   // parcel_cloud.density[passed_child_parcel_idx] = parcel_cloud.density[passed_parent_parcel_idx];
   // parcel_cloud.radius[passed_child_parcel_idx] = parcel_cloud.radius[passed_parent_parcel_idx];
   // parcel_cloud.density_tm1[passed_child_parcel_idx] = parcel_cloud.density_tm1[passed_parent_parcel_idx];

   // parcel_semi_mass_new = parcel_cloud.density[passed_child_parcel_idx] * parcel_cloud.radius[passed_child_parcel_idx] *
   //                        parcel_cloud.radius[passed_child_parcel_idx] * parcel_cloud.radius[passed_child_parcel_idx];

   // parcel_cloud.num_drop[passed_child_parcel_idx] =
   //     parcel_cloud.num_drop[passed_child_parcel_idx] * parcel_semi_mass_old / parcel_semi_mass_new;

       parcel_cloud.parcel_index[passed_child_parcel_idx] = parcel_cloud.parcel_index[passed_parent_parcel_idx];
       parcel_cloud.cloud_index[passed_child_parcel_idx] = parcel_cloud.cloud_index[passed_parent_parcel_idx];


 
   // printf("PARCEL_PROP.C L70 P_sat = %f\n",P_sat);
   // printf("PARCEL_PROP.C L67\n");

   // printf("\n ambient_pres = %f \n",ambient_pres);
   // printf("\n r_bubble old = %f \n", parcel_cloud.r_bubble[passed_parcel_idx]);


   //Need to reset properties but not flags for breakup (Assumes breakup starts again if it hasn't already occured)
   // parcel_cloud.r_bubble[passed_child_parcel_idx] = 2.0 * parcel_cloud.surf_ten[passed_parent_parcel_idx] / (P_sat - ambient_pres);
   if (parcel_cloud.r_bubble[passed_child_parcel_idx] < 0.0)
   {
      parcel_cloud.r_bubble[passed_child_parcel_idx] = 0.0;
      parcel_cloud.r_bubble_0[passed_child_parcel_idx] = 0.0;
   }
   parcel_cloud.r_bubble[passed_child_parcel_idx] = 0.0;

   // printf("\n PARCEL_PROP.C L69 r_bubble = %e\n", parcel_cloud.r_bubble[passed_parent_parcel_idx]);
   parcel_cloud.v_bubble[passed_child_parcel_idx] = 0.0;
   parcel_cloud.r_bubble_0[passed_child_parcel_idx] = 0.0;
   //Set these to prevent secondary thermal breakup 
   parcel_cloud.thermal_breakup_flag[passed_child_parcel_idx] = 5;
   parcel_cloud.pbt[passed_child_parcel_idx] = 0;
   // printf("\n END OF PARCEL_PROP.C \n");
   // printf("\n\n r_bubble = %e 	r_bubble_0 = %e", parcel_cloud.r_bubble[passed_parent_parcel_idx], parcel_cloud.r_bubble_0[passed_parent_parcel_idx]);

   // R_D_0
   parcel_cloud.r_drop_0[passed_child_parcel_idx] = parcel_cloud.radius[passed_parent_parcel_idx];

   // Droplet radius for thermal breakup purposes
   parcel_cloud.r_therm[passed_child_parcel_idx] = parcel_cloud.radius[passed_parent_parcel_idx];
   // Zero Omega and Eta (breakup variables)
   parcel_cloud.int_omega[passed_child_parcel_idx] = 0;
   parcel_cloud.omega[passed_child_parcel_idx] = parcel_cloud.omega_tm1[passed_parent_parcel_idx] = 0;
   parcel_cloud.eta_drop[passed_child_parcel_idx] = 0;
   parcel_cloud.dgre_cycle_count[passed_child_parcel_idx] = 0;

   // If parent's thermal_breakup_flag is set, displace the child parcel
   if (parcel_cloud.thermal_breakup_flag[passed_parent_parcel_idx] > 0)
   {
      // Get the velocity difference between child and parent
      CONVERGE_vec3_t velocity_diff, displacement;
      CONVERGE_vec3_diff(parcel_cloud.child_uu[passed_parent_parcel_idx], 
                       parcel_cloud.uu[passed_parent_parcel_idx], 
                       &velocity_diff);

      // Normalize velocity difference to get direction
      CONVERGE_vec3_normalize(velocity_diff);
      
      // Scale velocity_diff by parent's radius to get displacement vector
      CONVERGE_vec3_scale(velocity_diff, parcel_cloud.radius[passed_parent_parcel_idx]);
      
      // Use velocity_diff as displacement vector since it's already scaled

      CONVERGE_vec3_dup(velocity_diff, displacement);
      
      // Add displacement to parent's position to get child's position
      CONVERGE_vec3_add(parcel_cloud.xx[passed_parent_parcel_idx], displacement, &parcel_cloud.xx[passed_child_parcel_idx]);
      
      // Debug print (can be removed later)
      printf("\nParcel_Prop.c - Displacing parcel \n"
             "passed_parent_parcel_idx = %i\n"
             "passed_child_parcel_idx = %i\n"
             "parent_radius = %e\n"
             "velocity_diff = %e %e %e\n"
             "child_uu = %e %e %e\n"
             "child_uu address = %p\n"
             "displacement = %e %e %e\n"
             "child position = %e %e %e\n",
             passed_parent_parcel_idx,
             passed_child_parcel_idx,
             parcel_cloud.radius[passed_parent_parcel_idx],
             velocity_diff[0], velocity_diff[1], velocity_diff[2],
             parcel_cloud.child_uu[passed_parent_parcel_idx][0], parcel_cloud.child_uu[passed_parent_parcel_idx][1], parcel_cloud.child_uu[passed_parent_parcel_idx][2],
             (void*)parcel_cloud.child_uu[passed_parent_parcel_idx],
             displacement[0], displacement[1], displacement[2],
             parcel_cloud.xx[passed_child_parcel_idx][0], 
             parcel_cloud.xx[passed_child_parcel_idx][1], 
             parcel_cloud.xx[passed_child_parcel_idx][2]);
   }
   // }
      // parcel_cloud.tbt[passed_child_parcel_idx] = 0;


   parcel_cloud.m0[passed_child_parcel_idx] = (1.33333 * PI * CONVERGE_cube(parcel_cloud.radius[passed_parent_parcel_idx]) * parcel_cloud.num_drop[passed_parent_parcel_idx]);
   //Don't reset breakup flag
   // // printf("\n\n parcel_child called for parcel %i with parent %i which has temp %f\nradius %e and num_drop %e\nparent_radius = %e, parent_num_drop = %e\n",
   //        passed_child_parcel_idx, passed_parent_parcel_idx, parcel_cloud.temp[passed_child_parcel_idx],
   //        parcel_cloud.radius[passed_child_parcel_idx], parcel_cloud.num_drop[passed_child_parcel_idx],
   //        parcel_cloud.radius[passed_parent_parcel_idx], parcel_cloud.num_drop[passed_parent_parcel_idx]);
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
