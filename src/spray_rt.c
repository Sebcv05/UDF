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

/**********************************************************************/
/*                                                                    */
/* Name: spray_rt                                                     */
/*                                                                    */
/* Description: spray_rt performs the Rayleigh-Taylor (RT) drop       */
/* breakup calculation and includes viscosity effects on the breakup  */
/* time and size.                                                     */
/*                                                                    */
/* Variables passed in:  *parcel_cell   pointer to parcel cell        */
/*                       *parcel        pointer to spray parcel       */
/*                       *irt_break     pointer to flag indicating    */
/*                                      if RT breakup has taken place */
/*                                                                    */
/* Subroutines called: none                                           */
/*                                                                    */
/**********************************************************************/
static void init_tables(CONVERGE_species_t species);
static void destroy_tables(CONVERGE_species_t species);
static CONVERGE_table_t *pvap_table = NULL;

CONVERGE_UDF(spray_rt,
             IN(VALUE(CONVERGE_index_t, passed_parcel_idx),
                VALUE(CONVERGE_cloud_t, passed_spray_cloud),
                VALUE(CONVERGE_species_t, passed_species),
                VALUE(CONVERGE_precision_t, passed_density),
                VALUE(CONVERGE_precision_t, passed_mol_viscosity),
                VALUE(CONVERGE_precision_t *, passed_species_mass_fraction)),
             OUT(REF(CONVERGE_int_t, passed_irt_break)))

{
   struct ParcelCloud parcel_cloud;
   load_user_cloud(&parcel_cloud, passed_spray_cloud);

   init_tables(passed_species);

   CONVERGE_int_t iter, iter_max;
   CONVERGE_int_t rt_distribution_flag;
   CONVERGE_precision_t cnst3rt_parcel, cnst2b_parcel;
   CONVERGE_precision_t rt_accel, omax, kmax, wave_length;
   CONVERGE_precision_t breakup_time, rey_num;
   CONVERGE_precision_t term1, term2, term3, term4, term4_squared;
   CONVERGE_precision_t k_crit, k_low, k_high, dk1, dk;
   CONVERGE_precision_t radius_old, cd, converged, epsilon;
   CONVERGE_precision_t kk, kdk, omega_k, omega_kdk, domega;
   CONVERGE_precision_t delta_density, drop_smr, drop_radius;
   CONVERGE_precision_t random_num;
   CONVERGE_precision_t liquid_mass_in_cell;

   CONVERGE_index_t spray_elsa_flag   = CONVERGE_get_int("lagrangian.elsa_flag");
   CONVERGE_index_t num_gas_species   = CONVERGE_species_num_gas(passed_species);
   CONVERGE_index_t num_fluid_species = CONVERGE_species_num_fluid(passed_species);
   CONVERGE_precision_t dt            = CONVERGE_simulation_dt();

   // set RT model constants cnst3rt_parcel and cnst2b_parcel based on the injector
   // that the parcel came from

   CONVERGE_int_t from_injector = parcel_cloud.from_injector[passed_parcel_idx];

   //cnst3rt_parcel       = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_RT_CONST3);
   //cnst2b_parcel        = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_RT_CONST2);
   //rt_distribution_flag = CONVERGE_injector_get_parameter_flag(injector, INJECTOR_RT_DISTRIBUTION_FLAG);

   cnst3rt_parcel       = CONVERGE_spray_cloud_get_parameter_precision(passed_spray_cloud,PARCEL_IN_RT_SIZE_CONST);
   cnst2b_parcel        = CONVERGE_spray_cloud_get_parameter_precision(passed_spray_cloud,PARCEL_IN_RT_TIME_CONST);
   rt_distribution_flag = CONVERGE_spray_cloud_get_parameter_flag(passed_spray_cloud,PARCEL_IN_RT_DISTRIBUTION_FLAG);

   // initialize flag that determines if RT breakup has taken place to zero

   *passed_irt_break = 0;

   // Check for cooldown period after thermal breakup
   // If parcel underwent thermal breakup recently (lifetime < 10 microseconds), skip RT breakup
   CONVERGE_precision_t cooldown_time = 100.0e-6;  // 10 microseconds
   if(parcel_cloud.lifetime[passed_parcel_idx] < cooldown_time)  ////No RT breakup if parcel underwent thermal breakup 
   {
      destroy_tables(passed_species);
      return;
   }

   // calculate liquid viscosity and liquid reynolds number

   rey_num = 1.0e-10 + (2.0 * (parcel_cloud.radius[passed_parcel_idx]) * passed_density *
                        parcel_cloud.rel_vel_mag[passed_parcel_idx]) /
                          passed_mol_viscosity;

   // calculate drag coefficient including distortion effects

   if(rey_num > 1000.0)
   {
      cd = 0.424;
   }
   else
   {
      cd = (24.0 / rey_num) * (1.0 + (1.0 / 6.0) * CONVERGE_pow_2over3(rey_num));
   }

   cd = cd * (1.0 + 2.632 * parcel_cloud.distort[passed_parcel_idx]);

   // calculate acceleration

   rt_accel = (3.0 / 8.0) * cd * passed_density * parcel_cloud.rel_vel_mag[passed_parcel_idx] *
              parcel_cloud.rel_vel_mag[passed_parcel_idx] /
              (parcel_cloud.density[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx]);

   // perform bisection method search to determine the maximum growth rate and corresponding
   // wavenumber

   liquid_mass_in_cell = 0.0;
   if(spray_elsa_flag)
   {
      for(int isp = num_gas_species; isp < num_fluid_species; isp++)
      {
         liquid_mass_in_cell += passed_species_mass_fraction[isp] * passed_density;
      }
   }

   delta_density = fabs(parcel_cloud.density[passed_parcel_idx] - (passed_density - liquid_mass_in_cell)) + 1.0e-10;
   k_crit        = sqrt((delta_density)*rt_accel / parcel_cloud.surf_ten[passed_parcel_idx]);
   k_low         = 1.0e-10 * k_crit;
   k_high        = (0.95 * k_crit > 1.0e-99) ? 0.95 * k_crit : 1.0e-99;

   iter     = 0;
   iter_max = 500;

   dk1 = k_high * 1.0e-7;
   dk  = (1.0 < dk1) ? 1.0 : dk1;

   epsilon = 1.0e-4;

   converged = k_high - k_low;

   // calculate terms to go into growth rate expression

   term1 = (parcel_cloud.viscosity[passed_parcel_idx] + passed_mol_viscosity) /
           (parcel_cloud.density[passed_parcel_idx] + passed_density);
   term2 = rt_accel * (delta_density) / (parcel_cloud.density[passed_parcel_idx] + passed_density);
   term3 = parcel_cloud.surf_ten[passed_parcel_idx] / (parcel_cloud.density[passed_parcel_idx] + passed_density);
   term4 = (parcel_cloud.viscosity[passed_parcel_idx] + passed_mol_viscosity) /
           (parcel_cloud.density[passed_parcel_idx] + passed_density);
   term4_squared = term4 * term4;

   // iterate until convergence

   kk = 1.0e-10;
   while((iter < iter_max) && (converged >= epsilon))
   {
      iter = iter + 1;
      kk   = 0.5 * (k_low + k_high);

      kdk     = kk + dk;
      omega_k = -kk * kk * term1 + CONVERGE_sqrt(kk * term2 - kk * kk * kk * term3 + kk * kk * kk * kk * term4_squared);
      omega_kdk = -kdk * kdk * term1 +
                  CONVERGE_sqrt(kdk * term2 - kdk * kdk * kdk * term3 + kdk * kdk * kdk * kdk * term4_squared);

      domega = omega_kdk - omega_k;

      if(domega > 0.0)
      {
         k_low = kk;
      }
      else
      {
         k_high = kk;
      }

      converged = (k_high - k_low) / kk;
   }

   // calculate maximum growth rate and corresponding wavenumber
   //
   omax = -kk * kk * term1 + CONVERGE_sqrt(kk * term2 - kk * kk * kk * term3 + kk * kk * kk * kk * term4_squared);
   kmax = kk;

   // calculate wavelength corresponding to maximum growth rate

   wave_length = cnst3rt_parcel * 2.0 * PI / kmax;

   // RT waves are growing on the drop surface if the wavelegnth is smaller than the drop diameter

   if(wave_length < 2.0 * parcel_cloud.radius[passed_parcel_idx])
   {
      // calculate the breakup time and increment the time since last breakup

      breakup_time                              = cnst2b_parcel / omax;
      parcel_cloud.tbreak_rt[passed_parcel_idx] = parcel_cloud.tbreak_rt[passed_parcel_idx] + dt;

      if(parcel_cloud.tbreak_rt[passed_parcel_idx] > breakup_time)
      {
         // RT breakup takes place - update parcel radius and drop number
         CONVERGE_LAGRANGE_Spray_Parcel_make_small(passed_parcel_idx, passed_spray_cloud);
         radius_old  = parcel_cloud.radius[passed_parcel_idx];
         drop_smr    = 0.5 * wave_length;
         drop_radius = drop_smr;

         random_num  = CONVERGE_random_precision();
         drop_radius = CONVERGE_set_parcel_radius(rt_distribution_flag, drop_smr, random_num, from_injector);
         parcel_cloud.radius[passed_parcel_idx]     = drop_radius;
         parcel_cloud.radius_tm1[passed_parcel_idx] = parcel_cloud.radius[passed_parcel_idx];
         parcel_cloud.num_drop[passed_parcel_idx] =
            parcel_cloud.num_drop[passed_parcel_idx] * radius_old * radius_old * radius_old /
            (parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] *
             parcel_cloud.radius[passed_parcel_idx]);

         // since breakup has taken place, zero some of the parcel properties

         parcel_cloud.tbreak_rt[passed_parcel_idx]      = 0.0;
         parcel_cloud.distort[passed_parcel_idx]        = 0.0;
         parcel_cloud.distort_dot[passed_parcel_idx]    = 0.0;
         parcel_cloud.shed_num_drop[passed_parcel_idx]  = 0.0;
         parcel_cloud.shed_mass[passed_parcel_idx]      = 0.0;
         parcel_cloud.tbreak_kh[passed_parcel_idx]      = 0.0;
         parcel_cloud.tke0[passed_parcel_idx]           = 0.0;
         parcel_cloud.eps0[passed_parcel_idx]           = 0.0;
         parcel_cloud.lifetime[passed_parcel_idx]       = 0.0;
         parcel_cloud.area_reduction[passed_parcel_idx] = 0.0;

         // set flag that determines if RT breakup has taken place to 1

         *passed_irt_break = 1;
      }
   }

   // RT waves are not growing on the drop surface so zero the time since last breakup

   else
   {
      parcel_cloud.tbreak_rt[passed_parcel_idx] = 0.0;
   }

   destroy_tables(passed_species);

   return;
}


void init_tables(CONVERGE_species_t species)
{
   CONVERGE_iterator_t parcel_species_it;
   CONVERGE_species_parcel_iterator_create(species, &parcel_species_it);

   // Get the parcel species tables
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_PVAP_ID, &pvap_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
}

static void destroy_tables(CONVERGE_species_t species)
{
   CONVERGE_iterator_t parcel_species_it;
   CONVERGE_species_parcel_iterator_create(species, &parcel_species_it);

   // Destroy local parcel tables
   unload_species_tables(parcel_species_it, &pvap_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
}
