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
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Profiling variables
static int evap_call_count = 0;
static CONVERGE_precision_t total_evap_time = 0.0;
static CONVERGE_precision_t init_time = 0.0;
static CONVERGE_precision_t boil_time = 0.0;
static CONVERGE_precision_t evap_time = 0.0;
static CONVERGE_precision_t source_time = 0.0;
static CONVERGE_precision_t other_time = 0.0;  // Track unaccounted time

// Function to print profiling information
static void print_evap_profiling() {
 
    


    
    // Calculate unaccounted time
    CONVERGE_precision_t measured_time = init_time + boil_time + evap_time + source_time;
    other_time = (total_evap_time > measured_time) ? (total_evap_time - measured_time) : 0.0;
    
    printf("\n=== Spray Evaporation Profiling ===\n");
    printf("Total calls: %d\n", evap_call_count);
    printf("Total time: %.6f s (avg: %.6f ms/call)\n", 
           total_evap_time, (total_evap_time/(evap_call_count > 0 ? evap_call_count : 1))*1000.0);
    printf("Measured time: %.6f s (%.1f%% of total)\n", 
           measured_time, (measured_time/total_evap_time)*100.0);
    
    printf("\nTime distribution (of measured time):\n");
    printf("  Initialization: %8.2f%% (avg: %9.6f ms/call)\n", 
           100.0 * init_time / measured_time, (init_time/evap_call_count)*1000.0);
    printf("  Boiling calcs:  %8.2f%% (avg: %9.6f ms/call)\n", 
           100.0 * boil_time / measured_time, (boil_time/evap_call_count)*1000.0);
    printf("  Evaporation:    %8.2f%% (avg: %9.6f ms/call)\n", 
           100.0 * evap_time / measured_time, (evap_time/evap_call_count)*1000.0);
    printf("  Source terms:   %8.2f%% (avg: %9.6f ms/call)\n", 
           100.0 * source_time / measured_time, (source_time/evap_call_count)*1000.0);
    
    if (other_time > 0.0) {
        printf("  Unaccounted:    %8.2f%% (avg: %9.6f ms/call)\n", 
               100.0 * other_time / total_evap_time, (other_time/evap_call_count)*1000.0);
    }
    printf("=== End Profiling ===\n\n");
    fflush(stdout);  // Ensure output is flushed immediately
}

static void spray_evap_cell(CONVERGE_cloud_t cloud);
static void reset_parcel_temp_mfrac(const CONVERGE_index_t parcel_idx,
                                    CONVERGE_precision_t *temp_drop,
                                    CONVERGE_precision_t *temp_drop_starm1);
static CONVERGE_index_t compute_parcel_mole_fractions(CONVERGE_precision_t *mass_fraction,
                                                      CONVERGE_precision_t *mole_fraction,
                                                      CONVERGE_index_t num_parcel_species,
                                                      CONVERGE_species_t sp,
                                                      CONVERGE_index_t spray_evap_flag);


/********************************************************************************************/
/*                                                                                          */
/* Name: spray_evap,spray_evap_cell                                                         */
/*                                                                                          */
/* Description: spray_evap calculates spray droplet evaporation using an implicit scheme.   */
/*              spray_evap is a wrapper for spray_evap_cell.                                */
/*                                                                                          */
/* Variables passed in: none                                                                */
/*                                                                                          */
/* Subroutines called: spray_evap_cell                                                      */
/*                                                                                          */
/********************************************************************************************/

// Local UDF Field Variable Arrays
// Declared here to simplify function signature for spray_evap_cell
static const CONVERGE_precision_t *global_density;
static const CONVERGE_precision_t *global_mol_cond;
static const CONVERGE_precision_t *global_csubp;
static const CONVERGE_precision_t *global_csubv;
static const CONVERGE_precision_t *global_mol_viscosity;
static const CONVERGE_precision_t *global_temperature;
static const CONVERGE_precision_t *global_volume;
static const CONVERGE_precision_t *global_pressure;
static const CONVERGE_precision_t *global_sensible_sie;
static const CONVERGE_precision_t *global_species_massfrac;
static const CONVERGE_precision_t *global_species_massfrac_tm1;
static const CONVERGE_precision_t *global_density_tm1;
static CONVERGE_precision_t *global_src_ex_density;
static CONVERGE_precision_t *global_src_unburned_enth_evap_ex;
static CONVERGE_precision_t *global_src_sie_ex;
static CONVERGE_precision_t *global_src_species_ex;
static CONVERGE_precision_t *global_src_rif_fuel_evap_ex;
static CONVERGE_precision_t *global_ecfm_temp_evap_ex;
static const short *global_region_index;

static void init_tables(CONVERGE_species_t species);
static void destroy_tables(CONVERGE_species_t species);
static CONVERGE_table_t *pvap_table = NULL;
static CONVERGE_table_t *hvap_table = NULL;
static CONVERGE_table_t *visc_table = NULL;
static CONVERGE_table_t *cond_table = NULL;
static CONVERGE_table_t *cp_table   = NULL;
static CONVERGE_table_t *rho_table  = NULL;
static CONVERGE_table_t *h_table    = NULL;
static CONVERGE_table_t *evap_species_h_table;
static CONVERGE_table_t *evap_species_sensible_h_table;

static CONVERGE_species_t sp;
static CONVERGE_index_t num_gas_species;
static CONVERGE_index_t num_parcel_species;
static CONVERGE_index_t num_total_species;

static CONVERGE_precision_t dt;

static struct ParcelCloud parcel_cloud;

CONVERGE_UDF(spray_evap,
             IN(
                /*VELOCITY,*/
                FIELD(CONVERGE_precision_t *, density),
                FIELD(CONVERGE_precision_t *, mol_conductivity),
                VALUE(CONVERGE_mesh_t, mesh),
                FIELD(CONVERGE_vec3_t *, velocity),
                FIELD(CONVERGE_precision_t *, cp),
                FIELD(CONVERGE_precision_t *, cv),
                /*RAD_WEIGHT,*/
                /*SRC_RIF_FUEL_EVAP_EX*/
                /*VISCOSITY,*/
                FIELD(CONVERGE_precision_t *, mol_viscosity),
                FIELD(CONVERGE_precision_t *, temperature),
                FIELD(CONVERGE_precision_t *, volume),
                FIELD(CONVERGE_precision_t *, pressure),
                /*ENTHALPY,*/
                FIELD(CONVERGE_precision_t *, sensible_sie),
                ND_SCALAR_FIELD(species_mass_fraction),
                FIELD(CONVERGE_precision_t *, tm1_density),
                FIELD(CONVERGE_precision_t *, tm1_species_mass_fraction),
                /*SRC_EX_ENTHALPY,*/
                FIELD(short *, region_index)),
             OUT(FIELD(CONVERGE_precision_t *, src_ex_sensible_sie),
                 FIELD(CONVERGE_precision_t *, src_ex_species_mass_fraction),
                 FIELD(CONVERGE_precision_t *, src_ex_density),
                 FIELD(CONVERGE_precision_t *, ecfm3z_src_unburned_enth_evap_ex),
                 FIELD(CONVERGE_precision_t *, src_rif_fuel_evap_ex),
                 FIELD(CONVERGE_precision_t *, ecfm_temp_evap_ex)))
{
   CONVERGE_cloud_list_t spray_cloud_list = CONVERGE_mesh_get_spray_cloud_list(mesh);

   CONVERGE_iterator_t scl_it;
   CONVERGE_cloud_list_iterator_create(spray_cloud_list, &scl_it);
   for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(scl_it); i_pc != -1; i_pc = CONVERGE_iterator_next(scl_it))
   {
      CONVERGE_cloud_t cloud             = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, i_pc);
      CONVERGE_precision_t *src_evap_mom = CONVERGE_cloud_get_src_evap_mom(cloud);
      for(int jj = 0; jj < 3; jj++)
      {
         src_evap_mom[jj] = 0.0;
      }
   }

   CONVERGE_index_t spray_drag_flag = CONVERGE_get_int("lagrangian.drag_flag");

   // calculate relative velocity here if spray/gas coupling is turned off
   if(spray_drag_flag == 0)
   {
      CONVERGE_iterator_t pc_it;
      CONVERGE_cloud_t cloud;
      for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(scl_it); i_pc != -1; i_pc = CONVERGE_iterator_next(scl_it))
      {
         cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, i_pc);

         const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);
         const CONVERGE_vec3_t *uprime =
            (const CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(cloud, LAGRANGIAN_UPRIME);
         const CONVERGE_vec3_t *uu = (const CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(cloud, LAGRANGIAN_UU);
         CONVERGE_vec3_t *rel_vel  = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(cloud, LAGRANGIAN_REL_VEL);
         CONVERGE_precision_t *rel_vel_mag =
            (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cloud, LAGRANGIAN_REL_VEL_MAG);

         // zero out the evap source terms

         CONVERGE_cloud_iterator_create(cloud, &pc_it);
         for(CONVERGE_index_t i_p = CONVERGE_iterator_first(pc_it); i_p != -1; i_p = CONVERGE_iterator_next(pc_it))
         {
            rel_vel_mag[i_p] = 0.0;
            for(int ii = 0; ii < 3; ii++)
            {
               rel_vel[i_p][ii] = velocity[node_index][ii] + uprime[i_p][ii] - uu[i_p][ii];
            }
            rel_vel_mag[i_p] = CONVERGE_vec3_magnitude(rel_vel[i_p]);
         }

         CONVERGE_cloud_iterator_destroy(&pc_it);
      }
   }

   // Get high-level API containers
   sp                 = CONVERGE_mesh_species(mesh);
   num_gas_species    = CONVERGE_species_num_gas(sp);
   num_parcel_species = CONVERGE_species_num_parcel(sp);
   num_total_species  = species_mass_fraction_dimension;

   /*
    * Initialize variables and tables required for spray_evap_cell().
    */
   dt = CONVERGE_simulation_dt();

   init_tables(sp);

   // Note: These global variables are global to this file.
   global_density                   = density;
   global_mol_cond                  = mol_conductivity;
   global_csubp                     = cp;
   global_csubv                     = cv;
   global_mol_viscosity             = mol_viscosity;
   global_temperature               = temperature;
   global_volume                    = volume;
   global_pressure                  = pressure;
   global_sensible_sie              = sensible_sie;
   global_species_massfrac_tm1      = tm1_species_mass_fraction;
   global_species_massfrac          = species_mass_fraction;
   global_src_species_ex            = src_ex_species_mass_fraction;
   global_density_tm1               = tm1_density;
   global_src_ex_density            = src_ex_density;
   global_src_sie_ex                = src_ex_sensible_sie;
   global_src_unburned_enth_evap_ex = ecfm3z_src_unburned_enth_evap_ex;

   global_src_rif_fuel_evap_ex      = src_rif_fuel_evap_ex;
   global_ecfm_temp_evap_ex         = ecfm_temp_evap_ex;

   global_region_index              = region_index;

   // Get Simulation Meta data
   CONVERGE_precision_t dt = CONVERGE_simulation_dt();

   CONVERGE_cloud_t cloud;
   for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(scl_it); i_pc != -1; i_pc = CONVERGE_iterator_next(scl_it))
   {
      cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, i_pc);

      const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);
      CONVERGE_precision_t delmas       = global_src_ex_density[node_index] * dt / global_density[node_index];

      if(delmas > 0.2)
      {
         CONVERGE_logger_verbose("BEFORE spray_evap_cell DELMAS = %17.6e at ncyc= %ld", delmas, CONVERGE_ncyc());
      }

      spray_evap_cell(cloud);

      delmas = global_src_ex_density[node_index] * dt / global_density[node_index];
      if(delmas > 0.2)
      {
         CONVERGE_logger_verbose("AFTER spray_evap_cell DELMAS = %17.6e at ncyc= %ld", delmas, CONVERGE_ncyc());
      }
   }

   CONVERGE_iterator_destroy(&scl_it);

   destroy_tables(sp);
}

void spray_evap_cell(CONVERGE_cloud_t cloud)
{
   // Start timing the entire function
   CONVERGE_precision_t start_time = CONVERGE_mpi_wtime();
   CONVERGE_precision_t section_start, section_end;
   
   // Section 1: Initialization
   section_start = CONVERGE_mpi_wtime();
   
   // Setup parcel species counts and iterator
   CONVERGE_iterator_t psp_it;
   CONVERGE_species_parcel_iterator_create(sp, &psp_it);

   // Load user defined cloud wrapper
   load_user_cloud(&parcel_cloud, cloud);

   // Setup parcel cloud iterator for this cloud
   CONVERGE_iterator_t pc_it;
   CONVERGE_cloud_iterator_create(cloud, &pc_it);

   // Get relevent flags
   CONVERGE_flag_t parcel_boil_correlation_flag = CONVERGE_get_int("simulation.parcel_boil_correlation");
   CONVERGE_flag_t parcel_boil_flag;

   // Allocate local memory
   CONVERGE_index_t num_parcels        = CONVERGE_cloud_size(cloud);
   CONVERGE_precision_t *radius_new    = SAFE_calloc(num_parcels, CONVERGE_precision_t);
   CONVERGE_precision_t *moles         = SAFE_calloc(num_total_species, CONVERGE_precision_t);
   CONVERGE_precision_t *mol_frac      = SAFE_calloc(num_total_species, CONVERGE_precision_t);
   CONVERGE_precision_t *local_species = SAFE_calloc(num_total_species, CONVERGE_precision_t);

   CONVERGE_precision_t *vapor_mass_0               = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *vapor_mass                 = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *parcel_mole_fraction       = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *evap_min_radius            = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *evap_radius                = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *evap_mass_drop_0           = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *evap_mass_drop_1           = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *evap_cell_tot_evap_species = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *cell_tot_temp_species      = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   int *evap_all_flag                               = SAFE_calloc(num_parcel_species, int);
   int *parcel_species_boil_flag                    = SAFE_calloc(num_parcel_species, int);


   // Record initialization time
   section_end = CONVERGE_mpi_wtime();
   init_time += section_end - section_start;
   section_start = section_end;  // Reset section start for next timing section
   
   // Section 2: Boiling calculations
   section_start = CONVERGE_mpi_wtime();
   
   // Get spray_in variables


   // Get spray_in varables
   CONVERGE_index_t spray_evap_flag           = CONVERGE_get_int("lagrangian.evap_flag");
   CONVERGE_index_t drop_evap_source_flag     = CONVERGE_get_int("lagrangian.drop_evap_source_flag");
   CONVERGE_index_t urea_flag                 = CONVERGE_get_int("lagrangian.urea_flag");
   CONVERGE_index_t ecfm3z_flag               = CONVERGE_get_int("combust.ecfm_flag");
   
   CONVERGE_precision_t *region_energy_parcel_heating = CONVERGE_get_region_energy_parcel_heating();
   CONVERGE_precision_t *region_energy_parcel_phase_change = CONVERGE_get_region_energy_parcel_phase_change();
   CONVERGE_precision_t *region_energy_parcel_to_gas_phase = CONVERGE_get_region_energy_parcel_to_gas_phase();
   
   CONVERGE_precision_t spray_evap_n_diffuse  = CONVERGE_get_double("lagrangian.evap_n_diffuse");
   CONVERGE_precision_t spray_evap_d0_diffuse = CONVERGE_get_double("lagrangian.evap_d0_diffuse");
   CONVERGE_precision_t spray_scale_heat_trans_coeff_spray =
      CONVERGE_get_double("lagrangian.scale_heat_trans_coeff_spray");
   CONVERGE_precision_t spray_scale_mass_trans_coeff_spray =
      CONVERGE_get_double("lagrangian.scale_mass_trans_coeff_spray");
   
   CONVERGE_index_t hidden_multi_component_diffusion_flag = CONVERGE_get_int("hidden.multi_component_diffusion_flag");

   ///////////// Flash Boiling related local variables //////////////////////////
   CONVERGE_index_t evap_flag_flash_boiling = CONVERGE_get_int("lagrangian.flash_boiling.evap");
   CONVERGE_precision_t* volume_fraction = NULL;
   if( evap_flag_flash_boiling == 1 )
   {
      volume_fraction = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   }
   CONVERGE_precision_t distort_scale = CONVERGE_get_double("hidden.distort_scale");
   CONVERGE_precision_t evap_scale_factor_flash_boiling = CONVERGE_get_double("lagrangian.evap_scale_factor_flash_boiling");
   CONVERGE_precision_t pre_coeff_flshblg_l5 = CONVERGE_get_double("lagrangian.pre_coeff_flshblg_l5");
   CONVERGE_precision_t expnt_flshblg_l5 = CONVERGE_get_double("lagrangian.expnt_flshblg_l5");
   CONVERGE_precision_t pre_coeff_flshblg_l25 = CONVERGE_get_double("lagrangian.pre_coeff_flshblg_l25");
   CONVERGE_precision_t expnt_flshblg_l25 = CONVERGE_get_double("lagrangian.expnt_flshblg_l25");
   CONVERGE_precision_t pre_coeff_flshblg_g25 = CONVERGE_get_double("lagrangian.pre_coeff_flshblg_g25");
   CONVERGE_precision_t expnt_flshblg_g25 = CONVERGE_get_double("lagrangian.expnt_flshblg_g25");
   
   // Record boiling calculation time
   section_end = CONVERGE_mpi_wtime();
   boil_time += section_end - section_start;
   section_start = section_end;  // Reset section start for next timing section
   
   // Section 3: Evaporation calculations
   section_start = CONVERGE_mpi_wtime();
   // This section timing will be updated at the end of the section
   
   /////////////////////////////////////////////////////////////////////

   CONVERGE_precision_t *temp_boil = NULL;
   if( (parcel_boil_correlation_flag==1 && spray_evap_flag!=0) || (evap_flag_flash_boiling==1) )
   {
      temp_boil = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   }

   const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);

   const CONVERGE_precision_t min_spray_temp          = 200.0;
   const CONVERGE_precision_t min_spray_recovery_temp = 270.0;

   // Old table lookup vars
   CONVERGE_precision_t temp1;
   CONVERGE_precision_t temp2;
   CONVERGE_precision_t approx_sie;
   CONVERGE_precision_t average_hvap = 0.0;
   CONVERGE_precision_t cell_mass;
   CONVERGE_precision_t cell_tot_evap;
   CONVERGE_precision_t cond_term1 = 0.0;
   // CONVERGE_precision_t csubp_liquid1;
   CONVERGE_precision_t csubp_liquid;
   CONVERGE_precision_t delta_gas_temp;
   CONVERGE_precision_t delta_gas_temp_iterm1 = 0.0;
   CONVERGE_precision_t density1;
   CONVERGE_precision_t denth;
   CONVERGE_precision_t dmass;
   CONVERGE_precision_t drop_area;
   CONVERGE_precision_t drop_area_1;
   CONVERGE_precision_t enth;
   CONVERGE_precision_t enthalpy_liquid1;
   CONVERGE_precision_t enthalpy_liquid;
   CONVERGE_precision_t hvap;
   CONVERGE_precision_t local_density;
   CONVERGE_precision_t local_max_temp;
   CONVERGE_precision_t local_min_temp = 0.0;
   CONVERGE_precision_t local_sensible_sie;
   CONVERGE_precision_t local_sie;
   CONVERGE_precision_t local_source;
   CONVERGE_precision_t mass;
   CONVERGE_precision_t mass_0;
   CONVERGE_precision_t mass_drop_new;
   CONVERGE_precision_t mass_drop_tm1;
   CONVERGE_precision_t mass_star;
   CONVERGE_precision_t max_critical_temperature;
   CONVERGE_precision_t mol_visc = 0.0;
   CONVERGE_precision_t n_diffuse_m1;
   CONVERGE_precision_t new_temp_gas;
   CONVERGE_precision_t ro_mass_diff;
   CONVERGE_precision_t sc_num = 0.0;
   CONVERGE_precision_t src_sie_evap;
   CONVERGE_precision_t sum;
   CONVERGE_precision_t temp_gas_iterm1;
   CONVERGE_precision_t temp_new;
   CONVERGE_precision_t tg;
   CONVERGE_precision_t tot_evap_mass;
   CONVERGE_precision_t vapor_pres;
   CONVERGE_precision_t vaporization_term;
   CONVERGE_precision_t w_0, w_0_denom;

   CONVERGE_precision_t *mult_sc_num;
   CONVERGE_precision_t *mult_ro_mass_diff;
   CONVERGE_precision_t *mult_sh_num;
   //////////////////////Multi component diffusion local variables  ///////////////
   if(hidden_multi_component_diffusion_flag == 1)
   {
      mult_sc_num              = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
      mult_ro_mass_diff        = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
      mult_sh_num              = SAFE_calloc(num_parcels*num_parcel_species, CONVERGE_precision_t);
   }

   local_max_temp           = global_temperature[node_index];
   max_critical_temperature = 0.0;
   for(CONVERGE_index_t isp = CONVERGE_iterator_first(psp_it); isp != -1; isp = CONVERGE_iterator_next(psp_it))
   {
      const CONVERGE_precision_t isp_tcrit = CONVERGE_species_tcrit(sp, isp);   // CONVERGE_species_tcrit(sp,isp);
      if(isp_tcrit > local_max_temp)
      {
         local_max_temp = CONVERGE_species_tcrit(sp, isp);
      }
      if(isp_tcrit > max_critical_temperature)
      {
         max_critical_temperature = isp_tcrit;
      }
   }

   if( (parcel_boil_correlation_flag == 1 && spray_evap_flag != 0) || (evap_flag_flash_boiling==1) )
   {
      // boiling point (calculated based on the cell pressure)
      for(CONVERGE_index_t parcel_isp = 0, isp = CONVERGE_iterator_first(psp_it); isp != -1;
          parcel_isp++, isp                    = CONVERGE_iterator_next(psp_it))
      {
         const CONVERGE_precision_t tstep = 10.0; // Using a step of 10.0, this step is consistent with the input table
         const CONVERGE_precision_t isp_tcrit = CONVERGE_species_tcrit(sp, isp);
         CONVERGE_index_t boil_index          = (int)(isp_tcrit / tstep);
         CONVERGE_precision_t pvap1 = 0.0, pvap;
         while((pvap = CONVERGE_table_lookup(pvap_table[parcel_isp], boil_index * tstep)) > global_pressure[node_index])
         {
            boil_index--;
            pvap1 = pvap;
         }
         if(boil_index == (int)(isp_tcrit / tstep))
         {
            temp_boil[parcel_isp] = isp_tcrit;
         }
         else
         {
            temp_boil[parcel_isp] = tstep * (double)boil_index + tstep*(global_pressure[node_index] - pvap) / (pvap1 - pvap);
            temp_boil[parcel_isp] = fmin(temp_boil[parcel_isp], isp_tcrit);
         }
      }
   }

   local_max_temp = local_max_temp + 100.0;

   int local_evap_species_index = -1;
   if(drop_evap_source_flag == 0)
   {
      local_evap_species_index = CONVERGE_spray_evap_species_index();
   }

   // set initial values
   CONVERGE_precision_t temp_0 = global_temperature[node_index];
   
   // TODO: (rkratt) Change this in the internal code to make sense
   cell_mass = global_density[node_index] / (1.0 / global_volume[node_index]);   //FIXVOL

   // in the loop below we divide by the number of parcel species because we do not know the distribution of evaporating species in the gas phase.
   // this is only a problem for evap_source_flag=0.
   if(spray_evap_flag > 0)
   {
      if(drop_evap_source_flag == 0)
      {
         // ecfm3z model to be added
         for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
         {
            vapor_mass_0[isp] = global_species_massfrac[node_index * num_total_species + local_evap_species_index] *
                                cell_mass / ((CONVERGE_precision_t)num_parcel_species);
         }
      }
      else
      {
         for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
         {
            CONVERGE_index_t isp1 = CONVERGE_parcel_evap_species_lookup(isp);
            // ecfm3z model to be added
            // composite to be added
            vapor_mass_0[isp] = global_species_massfrac[node_index*num_total_species + isp1] * cell_mass;
         }
      }
   }

   memcpy(vapor_mass, vapor_mass_0, sizeof(CONVERGE_precision_t) * num_parcel_species);

   mass_0       = cell_mass;
   mass         = mass_0;
   n_diffuse_m1 = spray_evap_n_diffuse - 1.0;
   density1     = 1.293;

   CONVERGE_precision_t tot_mol = 0.0;
   for(int isp = 0; isp < num_gas_species; isp++)
   {
      moles[isp] =
         global_species_massfrac[node_index * num_total_species + isp] * cell_mass / CONVERGE_species_mw(sp, isp);
      tot_mol += moles[isp];
   }
   for(int isp = 0; isp < num_gas_species; isp++)
   {
      mol_frac[isp] = moles[isp] / tot_mol;
   }

   tot_mol = 0.0;
   w_0     = 0.0;
   w_0_denom = 1.0e-10;

   for(int isp = 0; isp < num_gas_species; isp++)
   {
      // ATTN: must skip gas species that are being sourced by parcel evaporation
      if(CONVERGE_skip_species_evap(isp) == 1)
      {
         continue;
      }
      w_0 += mol_frac[isp] * CONVERGE_species_mw(sp, isp);
      tot_mol += mol_frac[isp];
   }
   w_0 /= tot_mol;

   // implicit drop temperature solver
   CONVERGE_precision_t temp_gas = global_temperature[node_index];

   // calculate prandtl number followed by reynolds number, nusselt number and sherwood number
   CONVERGE_precision_t pr_num =
      global_mol_viscosity[node_index] * global_csubp[node_index] / global_mol_cond[node_index];

   // loop over all parcels in cell
   for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(pc_it); i_pc != -1; i_pc = CONVERGE_iterator_next(pc_it))
   {

      if((parcel_cloud.is_child[i_pc]==1 && parcel_cloud.lifetime[i_pc]<1.0e-4) || parcel_cloud.tbt[i_pc]){
         continue;
      }
      //printf("\n spray_evap_cell: L501, i_pc = %ld\n  ", i_pc);
      // see Borman and Ragland 1998 edition, p. 596
      tg       = (2.0 * parcel_cloud.temp[i_pc] + temp_gas) / 3.0;
      mol_visc = 0.0;
      for(int isp = 0; isp < num_gas_species; isp++)
      {
         mol_visc +=
            CONVERGE_table_lookup(visc_table[isp], tg) * global_species_massfrac[node_index * num_total_species + isp];
      }
      if(hidden_multi_component_diffusion_flag==1)
      {
         for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
         {
             mult_ro_mass_diff[isp] = density1*CONVERGE_get_D0_diffuse_species(isp)*pow((tg/273.0),(CONVERGE_get_n0_diffuse_species(isp)-1.0));
             mult_sc_num[isp] = mol_visc/mult_ro_mass_diff[isp];
         }
      }

      ro_mass_diff = density1 * spray_evap_d0_diffuse * pow((tg / 273.0), (n_diffuse_m1));
      sc_num       = mol_visc / ro_mass_diff;

      parcel_cloud.temp_starm1[i_pc]      = parcel_cloud.temp[i_pc];
      parcel_cloud.rey_num[i_pc]          = 2.0 * parcel_cloud.radius[i_pc] * global_density[node_index] *
                                      parcel_cloud.rel_vel_mag[i_pc] / (mol_visc + 1.0e-20) +
                                   1.0e-10;

      parcel_cloud.v_nu[i_pc] = 2.0 + 0.6 * sqrt(parcel_cloud.rey_num[i_pc]) * (CONVERGE_cbrt(pr_num));
      if(hidden_multi_component_diffusion_flag==1)
      {
         for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
         {
             mult_sh_num[i_pc * num_parcel_species + isp] = 2.0 + 0.6 * sqrt(parcel_cloud.rey_num[i_pc])*(CONVERGE_cbrt(mult_sc_num[isp]));
         }
      }
      parcel_cloud.v_sh[i_pc] = 2.0 + 0.6 * sqrt(parcel_cloud.rey_num[i_pc]) * (CONVERGE_cbrt(sc_num));
   }

   //******************************************************************************************************//
   //
   //       Start Implicit Solver (i.e., calculate drop temp and mass evaporated)
   //
   //******************************************************************************************************//
   const CONVERGE_precision_t local_volume = 1.0 / (1.0 / global_volume[node_index]);   //FIXVOL

   delta_gas_temp  = 1.0e9;
   temp_gas_iterm1 = 1.0e9;
   int iter_gas    = 0;
   while(fabs(delta_gas_temp) > 1.0 && fabs(temp_gas - temp_gas_iterm1) > 1.0)
   {
      iter_gas++;
      cell_tot_evap = 0.0;

      CONVERGE_precision_t bsub_d_avg = 0.0; /* initialized to prevent compiler warning */

      for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(pc_it); i_pc != -1; i_pc = CONVERGE_iterator_next(pc_it))
      {
         int inner_iter_flag                 = 1;
         int inner_iter                      = 0;
         CONVERGE_precision_t inner_iter_tol = 1.0;

         CONVERGE_precision_t tdrop        = parcel_cloud.temp[i_pc];
         CONVERGE_precision_t tdrop_starm1 = parcel_cloud.temp[i_pc];

         int max_inner_iter = 10;
         int min_inner_iter = 1;

         // omega is the under-relaxation. if convergence not working this can be lowered.
         CONVERGE_precision_t omega = 1.0;

         int recovery_flag    = 0;
         int recovery_counter = 0;

         if( evap_flag_flash_boiling == 1 )
         {
            double sum_volume_fraction = 0;
            for( CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++ )
            {
               double density_sp = CONVERGE_table_lookup(rho_table[isp], parcel_cloud.temp[i_pc]);

               volume_fraction[isp] = parcel_cloud.density[i_pc] * parcel_cloud.mfrac[i_pc*num_parcel_species+isp] / density_sp;

               sum_volume_fraction += volume_fraction[isp];
            }
            for( CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++ )
            {
               if( sum_volume_fraction < 1.e-27 )
               {
                  volume_fraction[isp] = 0.;
               }
               else
               {
                  volume_fraction[isp] = volume_fraction[isp] / sum_volume_fraction;
               }
            }
         }

         parcel_boil_flag = 0;
         if(parcel_boil_correlation_flag == 1 && spray_evap_flag != 0)
         {
            for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
            {
               if(temp_boil[isp] <= parcel_cloud.temp[i_pc] && temp_boil[isp] < temp_gas &&
                  parcel_cloud.temp[i_pc] < temp_gas &&
                  parcel_cloud.mfrac[i_pc * num_parcel_species + isp] >
                     0.0)   //0.01 is to indicate that the species has almost been removed
               {
                  parcel_boil_flag = 1;
                  parcel_species_boil_flag[isp] = 1;
               }
            }
         }

         while((inner_iter_flag == 1) && (inner_iter < max_inner_iter))
         {
            inner_iter_flag = 0;
            inner_iter += 1;

            int do_recovery = 0;

            // check if any component of the parcel is boiling. If so, temperature ofthe parcel will remain constant
            if( evap_flag_flash_boiling == 1 )
            {
               parcel_boil_flag = 0;
            }

            temp1 = tdrop;
            temp2 = (2.0 * tdrop + temp_gas) / 3.0;

            CONVERGE_precision_t thermal_cond_gas = 0.0;
            mol_visc                              = 0.0;
            for(int isp = 0; isp < num_gas_species; isp++)
            {
               thermal_cond_gas += global_species_massfrac[node_index * num_total_species + isp] *
                                   CONVERGE_table_lookup(cond_table[isp], temp2);

               mol_visc += global_species_massfrac[node_index * num_total_species + isp] *
                           CONVERGE_table_lookup(visc_table[isp], temp2);
            }

            ro_mass_diff = density1 * spray_evap_d0_diffuse * pow((temp2 / 273.0), (n_diffuse_m1));
            sc_num       = mol_visc / ro_mass_diff;

            if(hidden_multi_component_diffusion_flag==1)
            {
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                   mult_ro_mass_diff[isp] = density1*CONVERGE_get_D0_diffuse_species(isp)*pow((temp2/273.0),(CONVERGE_get_n0_diffuse_species(isp)-1.0));
                   mult_sc_num[isp] = mol_visc/mult_ro_mass_diff[isp];
               }
            }

            mass_drop_tm1 =
               parcel_cloud.density_tm1[i_pc] * (4.0 / 3.0) * PI * (CONVERGE_cube(parcel_cloud.radius_tm1[i_pc]));

            if((drop_evap_source_flag == 1 || drop_evap_source_flag == 2) && spray_evap_flag > 0)
            {
               compute_parcel_mole_fractions(&(parcel_cloud.mfrac[i_pc * num_parcel_species]), parcel_mole_fraction,
                       num_parcel_species, sp, spray_evap_flag);
               w_0_denom = w_0;
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  vapor_pres = CONVERGE_table_lookup(pvap_table[isp], temp1) + 1.0e-3;
                  vapor_pres = (vapor_pres > global_pressure[node_index]) ? (global_pressure[node_index]) : (vapor_pres);
                  // To be added: Urea model

                  w_0_denom += parcel_mole_fraction[isp] * vapor_pres *
                          CONVERGE_get_parcel_species_mw(isp)  / global_pressure[node_index] -
                          w_0 * parcel_mole_fraction[isp] * vapor_pres / global_pressure[node_index];
               }
            }

            if(parcel_boil_flag)
            {
               average_hvap = 0.0;
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  average_hvap += parcel_cloud.mfrac[i_pc * num_parcel_species + isp] *
                                  (CONVERGE_table_lookup(hvap_table[isp], temp1) + 1.0e-3);
               }
               if((1 + global_csubp[node_index] * (temp_gas - tdrop) / average_hvap) <= 0.0)
               {
                  parcel_boil_flag = 0;  //NaN avoider
               }
            }

            vaporization_term = 0.0;
            csubp_liquid      = 0.0;
            bsub_d_avg        = 0.0;
            for(CONVERGE_index_t isp = 0, isp1 = CONVERGE_iterator_first(psp_it); isp < num_parcel_species;
                isp++, isp1                    = CONVERGE_iterator_next(psp_it))
            {
               evap_all_flag[isp] = 0;
               if(parcel_cloud.mfrac_tm1[i_pc * num_parcel_species + isp] <= 0.0)
               {
                  parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = 0.0;
                  evap_mass_drop_0[isp]                               = 0.0;
                  continue;
               }

               evap_mass_drop_0[isp] = mass_drop_tm1 * parcel_cloud.mfrac_tm1[i_pc * num_parcel_species + isp];

               evap_min_radius[isp] =
                  (3.0 / 4.0) * (mass_drop_tm1 - evap_mass_drop_0[isp]) / (PI * parcel_cloud.density[i_pc]);
               evap_min_radius[isp] = (evap_min_radius[isp] > 1.0e-25) ? (CONVERGE_cbrt(evap_min_radius[isp])) : (0.0);

               vapor_pres = CONVERGE_table_lookup(pvap_table[isp], temp1) + 1.0e-3;
               vapor_pres = (vapor_pres > global_pressure[node_index]) ? (global_pressure[node_index]) : (vapor_pres);
               // Urea model to be added

               hvap = CONVERGE_table_lookup(hvap_table[isp], temp1) + 1.0e-3;

               csubp_liquid +=
                  parcel_cloud.mfrac[i_pc * num_parcel_species + isp] * CONVERGE_table_lookup(cp_table[isp], temp1);

               //skip the evap calc.need some variables above, such as csubp_liquid, to do the following calculation even if evap_flag==0
               if(spray_evap_flag == 0)
               {
                  parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = 0.0;
                  evap_mass_drop_0[isp]                               = 0.0;
                  continue;
               }

               CONVERGE_precision_t y1_star = 0.0;
               if(spray_evap_flag > 0)
               {
                  if(drop_evap_source_flag == 0)
                  {
                     y1_star = CONVERGE_species_mw(sp, local_evap_species_index) /
                               (CONVERGE_species_mw(sp, local_evap_species_index) +
                                (w_0 * (global_pressure[node_index] / vapor_pres - 1.0)));
                  }
                  else
                  {
                     // Urea model to be implemented
                     y1_star = parcel_mole_fraction[isp] * vapor_pres /
                             global_pressure[node_index] * CONVERGE_get_parcel_species_mw(isp) / w_0_denom;
                  }
               }
               y1_star                     = (y1_star < 0.0) ? 1.0e-10 : y1_star;
               y1_star                     = (y1_star < 1.0) ? y1_star : 1.0;
               CONVERGE_precision_t y1     = vapor_mass[isp] / mass;
               CONVERGE_precision_t bsub_d = (y1_star - y1) / (1.0 - y1_star + 1.0e-10);

               //   Only allow evap, No condensation
               bsub_d = (bsub_d < 1.0e-15) ? 1.0e-15 : bsub_d;
               bsub_d_avg += bsub_d * parcel_cloud.mfrac[i_pc * num_parcel_species + isp];
               CONVERGE_precision_t log_bsub_d = log(bsub_d + 1.0);


               //printf("\n spray_evap_cell: L750, bsub_d = %e\n  ", bsub_d);

               CONVERGE_precision_t mass_trans_coeff = 0.0;
               if(hidden_multi_component_diffusion_flag == 1)
               {
                   mass_trans_coeff = spray_scale_mass_trans_coeff_spray * mult_sh_num[i_pc * num_parcel_species + isp] * (mol_visc / mult_sc_num[isp]) /
                                   (2.0 * parcel_cloud.radius[i_pc] * parcel_cloud.density[i_pc]);
               }
               else
               {
                  mass_trans_coeff = spray_scale_mass_trans_coeff_spray * parcel_cloud.v_sh[i_pc] * (mol_visc / sc_num) /
                                  (2.0 * parcel_cloud.radius[i_pc] * parcel_cloud.density[i_pc]);
               }

               if(parcel_boil_flag == 1 && parcel_species_boil_flag[isp] == 1)   // boiling correlation
               {
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] =
                     -global_mol_cond[node_index] /
                     (parcel_cloud.density[i_pc] * global_csubp[node_index] * parcel_cloud.radius[i_pc]) *
                     (1 + 0.23 * sqrt(parcel_cloud.rey_num[i_pc])) *
                     log(1 + global_csubp[node_index] * (temp_gas - tdrop) / average_hvap);
               }
               else
               {
                  if(spray_evap_flag == 1)
                  {
                     parcel_cloud.drdt[i_pc * num_parcel_species + isp] = -mass_trans_coeff * log_bsub_d;
                  }

                  if(spray_evap_flag == 2)
                  {
                     parcel_cloud.drdt[i_pc * num_parcel_species + isp] =
                        -mass_trans_coeff * (y1_star - y1) * (bsub_d / (pow((1.0 + bsub_d), 0.568)));
                  }
               }
               //printf("\n spray_evap_cell: L785, parcel_cloud.drdt[i_pc * num_parcel_species + isp] = %e\n  ", parcel_cloud.drdt[i_pc * num_parcel_species + isp]);
               if( evap_flag_flash_boiling==1 )
               {
                  //printf("\n is child %d\n", parcel_cloud.is_child[i_pc]);
                  //printf("\n lifetime %e\n", parcel_cloud.lifetime[i_pc]);
                  //printf("\n tdrop %e\n", tdrop);
                  //printf("\n temp_boil %e\n", temp_boil[isp]);
                  double super_heat_degree = tdrop - temp_boil[isp];
                  //printf("\n spray_evap_cell: L789, super_heat_degree = %e\n  ", super_heat_degree);

                  if( super_heat_degree > 0.2 )
                  {
                     double density_sp = CONVERGE_table_lookup(rho_table[isp], tdrop);
                     if( super_heat_degree < 5 )
                     {
                        parcel_cloud.drdt[i_pc * num_parcel_species + isp] -= 
                           evap_scale_factor_flash_boiling * volume_fraction[isp] *
                           pre_coeff_flshblg_l5 * pow(super_heat_degree, expnt_flshblg_l5) *
                           super_heat_degree / density_sp / hvap;
                     }
                     else if( super_heat_degree < 25 && super_heat_degree >= 5 )
                     {
                        parcel_cloud.drdt[i_pc * num_parcel_species + isp] -= 
                           (1.0 + distort_scale*parcel_cloud.distort[i_pc]) *
                           evap_scale_factor_flash_boiling * volume_fraction[isp] *
                           pre_coeff_flshblg_l25 * pow(super_heat_degree, expnt_flshblg_l25) *
                           super_heat_degree / density_sp / hvap;
                     }
                     else
                     {
                        parcel_cloud.drdt[i_pc * num_parcel_species + isp] -=
                           (1.0 + distort_scale * parcel_cloud.distort[i_pc]) *
                           evap_scale_factor_flash_boiling * volume_fraction[isp] *
                           pre_coeff_flshblg_g25 * pow(super_heat_degree, expnt_flshblg_g25) *
                           super_heat_degree / density_sp / hvap;
                     }
                  }
               }
               //printf("\n spray_evap_cell L815 ");
               if(parcel_cloud.is_child[i_pc]==2)
               {
                  //printf("\n spray_evap_cell: parcel is child\n");
                  if(parcel_cloud.lifetime[i_pc] < 1.0e-4)
                  {
                     //printf("\n tdrop = %e\n ", tdrop);
                     //printf("\n temp_boil = %e\n", temp_boil[isp]);
                     double super_heat_degree = tdrop - temp_boil[isp];

                     //printf("\n spray_evap_cell: parcel lifetime < 1.0e-4\n");   
                     // parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
                     double density_sp = CONVERGE_table_lookup(rho_table[isp], tdrop);
                     if( super_heat_degree < 5 && super_heat_degree >0.2)
                     {
                        //printf("super_heat_degree < 5\n");
                        parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
                        // parcel_cloud.drdt[i_pc * num_parcel_species + isp] =- 
                        //    evap_scale_factor_flash_boiling * volume_fraction[isp] *
                        //    pre_coeff_flshblg_l5 * pow(super_heat_degree, expnt_flshblg_l5) *
                        //    super_heat_degree / density_sp / hvap;
                     }
                     else if( super_heat_degree < 25 && super_heat_degree >= 5 )
                     {
                        // printf("super_heat_degree < 25 && super_heat_degree >= 5\n");
                        parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
                        // parcel_cloud.drdt[i_pc * num_parcel_species + isp] =- 
                        //    (1.0 + distort_scale*parcel_cloud.distort[i_pc]) *
                        //    evap_scale_factor_flash_boiling * volume_fraction[isp] *
                        //    pre_coeff_flshblg_l25 * pow(super_heat_degree, expnt_flshblg_l25) *
                        //    super_heat_degree / density_sp / hvap;
                     }
                     else
                     {
                        //printf("super_heat_degree >= 25\n");
                        parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
                        // parcel_cloud.drdt[i_pc * num_parcel_species + isp] =-
                        //    (1.0 + distort_scale * parcel_cloud.distort[i_pc]) *
                        //    evap_scale_factor_flash_boiling * volume_fraction[isp] *
                        //    pre_coeff_flshblg_g25 * pow(super_heat_degree, expnt_flshblg_g25) *
                        //    super_heat_degree / density_sp / hvap;
                     }
                  }
            }

               if((parcel_cloud.drdt[i_pc * num_parcel_species + isp] * dt + parcel_cloud.radius[i_pc]) <
                  evap_min_radius[isp])
               {
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] =
                     -((parcel_cloud.radius[i_pc] - evap_min_radius[isp]) / dt);
                  evap_all_flag[isp] = 1;
               }

               //  again don't allow condensation
               if(parcel_cloud.drdt[i_pc * num_parcel_species + isp] >= 0.0)
               {
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
               }

               // temperature can't be above critical temp
               const CONVERGE_precision_t isp_tcrit = CONVERGE_species_tcrit(sp, isp1);
               if(tdrop >= (isp_tcrit - 1.0e-6))   // make the species disappear
               {
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] =
                     -((parcel_cloud.radius[i_pc] - evap_min_radius[isp]) / dt);
                  evap_all_flag[isp] = 1;
               }

               evap_radius[isp] = parcel_cloud.radius[i_pc] + (parcel_cloud.drdt[i_pc * num_parcel_species + isp] * dt);
               evap_mass_drop_1[isp] =
                  (4.0 / 3.0) * PI * parcel_cloud.density[i_pc] * (CONVERGE_cube(evap_radius[isp]));
               parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = (evap_mass_drop_1[isp] - mass_drop_tm1) / dt;

               // Do not evaporate more mass than available
               if((-parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt) > evap_mass_drop_0[isp])
               {
                  parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = -evap_mass_drop_0[isp] / dt;
               }
               vaporization_term =
                  vaporization_term + (parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt) * hvap;
            }

            mass_drop_new = mass_drop_tm1;
            for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
            {
               mass_drop_new += parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt;
            }

            if(mass_drop_new > 1.0e-36) /* normal update */
            {
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  parcel_cloud.mfrac[i_pc * num_parcel_species + isp] =
                     (evap_mass_drop_0[isp] + parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt) /
                     (mass_drop_new);
                  if(parcel_cloud.mfrac[i_pc * num_parcel_species + isp] < 0.0 || evap_all_flag[isp] == 1)
                  {
                     parcel_cloud.mfrac[i_pc * num_parcel_species + isp] = 0.0;
                  }
               }
            }
            else /* if mass drop less than 1.0e-36 use old mass fractions */
            {
               mass_drop_new = 1.0e-36;
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  parcel_cloud.mfrac[i_pc * num_parcel_species + isp] =
                     parcel_cloud.mfrac_tm1[i_pc * num_parcel_species + isp];
               }
            }

            sum = 0.0;
            for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
            {
               sum += parcel_cloud.mfrac[i_pc * num_parcel_species + isp];
            }
            if(sum > 0.1) /* normalize */
            {
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  parcel_cloud.mfrac[i_pc * num_parcel_species + isp] =
                     parcel_cloud.mfrac[i_pc * num_parcel_species + isp] / sum;
               }
            }
            else /* all species evaporated */
            {
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  parcel_cloud.mfrac[i_pc * num_parcel_species + isp] =
                     parcel_cloud.mfrac_tm1[i_pc * num_parcel_species + isp];
               }
            }

            radius_new[i_pc] = (3.0 / 4.0) * mass_drop_new / (PI * parcel_cloud.density[i_pc]);
            radius_new[i_pc] = CONVERGE_cbrt(radius_new[i_pc]);
            radius_new[i_pc] =
               fmin(radius_new[i_pc],
                    parcel_cloud.radius[i_pc]);   //Truncation errors for very small drdt make this necessary

            drop_area   = 4.0 * PI * parcel_cloud.radius[i_pc] * parcel_cloud.radius[i_pc];
            drop_area_1 = 4.0 * PI * radius_new[i_pc] * radius_new[i_pc];

            drop_area = 0.5 * (drop_area + drop_area_1);

            CONVERGE_precision_t spray_wall_heat_source = 0.0;   // Wruck spray-wall heat exchange

            CONVERGE_precision_t heat_trans_coeff;
            heat_trans_coeff = spray_scale_heat_trans_coeff_spray * parcel_cloud.v_nu[i_pc] * thermal_cond_gas /
                               (parcel_cloud.radius[i_pc] + radius_new[i_pc]);

            if(parcel_boil_flag)   // boiling
            {
               // this is to make tdrop same as parcel's temperature
               cond_term1 = (fabs(tdrop - temp_gas) == 0.0) ? (0.0) : ((vaporization_term) / (tdrop - temp_gas));
            }
            else
            {
               if(spray_evap_flag == 1)
               {
                  cond_term1 =
                     (bsub_d_avg == 0.0) ? 0.0 : dt * drop_area * heat_trans_coeff * log(1.0 + bsub_d_avg) / bsub_d_avg;
               }
               else if(spray_evap_flag == 2)
               {
                  cond_term1 = dt * drop_area_1 * heat_trans_coeff * (1.0 / (pow((1.0 + bsub_d_avg), 0.678)));
               }
               else if(spray_evap_flag == 0)
               {
                  cond_term1 = dt * drop_area * heat_trans_coeff;
               }
            }
            //Turn of spalding number correlation for children after breakup 
            if(parcel_cloud.is_child[i_pc])
            {
               if(parcel_cloud.lifetime[i_pc] < 1.0e-4)
               {
                  // cond_term1 = dt * drop_area * heat_trans_coeff;
                  cond_term1 = 0.0; //Zero
                  // if(i_pc%100==0){
                  // // printf(" setting cond_term1 to zero for child %ld\n", i_pc);
                  // }
               }
            }

            CONVERGE_precision_t denom = cond_term1 + (csubp_liquid * mass_drop_new);

            tdrop = (denom == 0.0) ? tdrop_starm1
                                   : ((csubp_liquid * mass_drop_new * parcel_cloud.temp[i_pc]) + vaporization_term +
                                      (cond_term1 * temp_gas) + spray_wall_heat_source) /
                                        denom;

            tdrop = omega * (tdrop - tdrop_starm1) + tdrop_starm1;

            if((tdrop < min_spray_recovery_temp))
            {
               tdrop       = (tdrop < min_spray_temp) ? (min_spray_temp) : (tdrop);
               do_recovery = 1;
            }
            if(tdrop >= max_critical_temperature)
            {
               tdrop = max_critical_temperature;
               if(temp_gas < max_critical_temperature)
               {
                  do_recovery = 1;
               }
            }

            if((parcel_boil_correlation_flag && spray_evap_flag != 0) && evap_flag_flash_boiling==0)
            {
               parcel_boil_flag = 0;
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  parcel_species_boil_flag[isp] = 0;
                  if(tdrop >= temp_boil[isp] && parcel_cloud.mfrac[i_pc * num_parcel_species + isp] > 0.0)
                  {
                     parcel_boil_flag = 1;
                     parcel_species_boil_flag[isp] = 1;
                     inner_iter_flag = 1;
                     min_inner_iter += 1;
                  }
               }
            }

            // This is to correct omega if new solution goes out of bounds.
            if(do_recovery)
            {
               recovery_flag   = 1;
               inner_iter_flag = 1;
               recovery_counter += 1;
               omega *= 0.5;
               inner_iter_tol *= 0.5;
               min_inner_iter *= 2;

               if(recovery_counter == 1)
               {
                  max_inner_iter = 40;
               }

               if(recovery_counter > 3)
               {
                  tdrop = tdrop_starm1;
                  break;
               }

               reset_parcel_temp_mfrac(i_pc, &tdrop, &tdrop_starm1);
            }

            // This checks if solution is converged.
            if(fabs(tdrop - tdrop_starm1) > inner_iter_tol && tdrop > min_spray_temp)
            {
               inner_iter_flag = 1;

               if((inner_iter == 10) && !recovery_flag)
               {
                  max_inner_iter = 20;
                  omega *= 0.5;
                  inner_iter_tol *= 0.5;
                  min_inner_iter *= 2;

                  reset_parcel_temp_mfrac(i_pc, &tdrop, &tdrop_starm1);
               }
            }

            inner_iter_flag = (inner_iter < min_inner_iter) ? (1) : (inner_iter_flag);

            tdrop_starm1 = tdrop;
         }

         parcel_cloud.temp_starm1[i_pc] = tdrop;
      }

      for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
      {
         vapor_mass[isp] = vapor_mass_0[isp];
      }
      mass = mass_0;

      for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(pc_it); i_pc != -1; i_pc = CONVERGE_iterator_next(pc_it))
      {
         dmass = 0.0;
         for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
         {
            dmass -= parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt * parcel_cloud.num_drop[i_pc];
            vapor_mass[isp] -= parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt * parcel_cloud.num_drop[i_pc];
         }

         mass = mass + dmass;
      }

      //    update gas temperature
      denth          = 0.0;
      local_min_temp = global_temperature[node_index];
      for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(pc_it); i_pc != -1; i_pc = CONVERGE_iterator_next(pc_it))
      {
         tot_evap_mass = 0.0;
         for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
         {
            tot_evap_mass =
               tot_evap_mass - parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt * parcel_cloud.num_drop[i_pc];
         }
         cell_tot_evap = cell_tot_evap + tot_evap_mass;

         {
            temp1 = parcel_cloud.temp_tm1[i_pc];
            temp2 = parcel_cloud.temp_starm1[i_pc];

            local_min_temp =
               (local_min_temp > parcel_cloud.temp_starm1[i_pc]) ? (parcel_cloud.temp_starm1[i_pc]) : (local_min_temp);

            csubp_liquid = 0.0;
            //csubp_liquid1    = 0.0;
            enthalpy_liquid  = 0.0;
            enthalpy_liquid1 = 0.0;
            hvap             = 0.0;

            for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
            {
               csubp_liquid +=
                  CONVERGE_table_lookup(cp_table[isp], temp1) * parcel_cloud.mfrac[i_pc * num_parcel_species + isp];
               //csubp_liquid1 =
               //   CONVERGE_table_lookup(cp_table[isp], temp2) * parcel_cloud.mfrac[i_pc * num_parcel_species + isp];
               enthalpy_liquid +=
                  CONVERGE_table_lookup(h_table[isp], temp1) * parcel_cloud.mfrac[i_pc * num_parcel_species + isp];
               enthalpy_liquid1 +=
                  CONVERGE_table_lookup(h_table[isp], temp2) * parcel_cloud.mfrac[i_pc * num_parcel_species + isp];
               hvap -= CONVERGE_table_lookup(hvap_table[isp], temp2) *
                       parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt * parcel_cloud.num_drop[i_pc];
            }

            enth = 0.0;
            if(spray_evap_flag > 0)
            {
               if(drop_evap_source_flag == 0)
               {
                  enth = CONVERGE_table_lookup(evap_species_sensible_h_table[local_evap_species_index], temp1) * tot_evap_mass;
               }
               else
               {
                  for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
                  {
                     // Urea model to be implemented
                     CONVERGE_index_t isp1 = CONVERGE_parcel_evap_species_lookup(isp);
                     enth += CONVERGE_table_lookup(evap_species_sensible_h_table[isp1], temp1) *
                             (-parcel_cloud.dm_dt[i_pc * num_parcel_species + isp]) * dt * parcel_cloud.num_drop[i_pc];
                  }
               }
            }

            //   radius can't be less than 0.0
            if(radius_new[i_pc] < 0.0)
            {
               radius_new[i_pc] = 0.0;
            }
            mass_star = (4.0 / 3.0) * PI * parcel_cloud.density[i_pc] * parcel_cloud.num_drop[i_pc] * radius_new[i_pc] *
                        radius_new[i_pc] * radius_new[i_pc];

            //  heat of vaporization should be added to the evaporated liquid
            denth += enth - hvap - mass_star * (enthalpy_liquid1 - enthalpy_liquid);
         }
      }

      // Approximate gas side temp this is where we make our next guess
      delta_gas_temp =
         (temp_0 + denth / ((local_volume * global_density[node_index] + cell_tot_evap) * global_csubp[node_index])) -
         temp_gas;

      if(iter_gas == 1)   // we don't have two points, yet, use the delta_gas_temp for temp_gas
      {
         new_temp_gas = temp_gas + delta_gas_temp;
      }
      else if(fabs(delta_gas_temp) < 1.0 || fabs(temp_gas - temp_gas_iterm1) < 1.0 ||
              fabs(delta_gas_temp - delta_gas_temp_iterm1) < 0.01)
      {
         // don't need to do anything, the loop is converged
         new_temp_gas = temp_gas;
      }
      else
      {
         // calculate the new temp_gas using Newton's method
         new_temp_gas =
            temp_gas + delta_gas_temp * ((temp_gas_iterm1 - temp_gas) / (delta_gas_temp - delta_gas_temp_iterm1));
      }

      temp_gas_iterm1       = temp_gas;
      delta_gas_temp_iterm1 = delta_gas_temp;

      temp_gas = new_temp_gas;

      if(temp_gas < local_min_temp)
      {
         temp_gas = local_min_temp;
      }
      if(temp_gas > local_max_temp)
      {
         temp_gas = local_max_temp;
      }

      if(iter_gas > 20)
      {
         break;
      }
   }
   // *********************************************************************************************************** //
   //
   //       End Implicit Solver
   //
   // *********************************************************************************************************** //

   // update parcel radius
   mass = local_volume * global_density[node_index];

   denth         = 0.0;
   cell_tot_evap = 0.0;
   memset(evap_cell_tot_evap_species, 0, sizeof(CONVERGE_precision_t) * num_parcel_species);
   memset(cell_tot_temp_species, 0, sizeof(CONVERGE_precision_t) * num_parcel_species);

   for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(pc_it); i_pc != -1; i_pc = CONVERGE_iterator_next(pc_it))
   {
      tot_evap_mass = 0.0;
      for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
      {
         tot_evap_mass -= parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt * parcel_cloud.num_drop[i_pc];
      }

      mass_star = (4.0 / 3.0) * PI * parcel_cloud.density[i_pc] * parcel_cloud.num_drop[i_pc] * radius_new[i_pc] *
                  radius_new[i_pc] * radius_new[i_pc];
      dmass = (4.0 / 3.0) * PI * parcel_cloud.density[i_pc] * parcel_cloud.num_drop[i_pc] *
                 parcel_cloud.radius_tm1[i_pc] * parcel_cloud.radius_tm1[i_pc] * parcel_cloud.radius_tm1[i_pc] -
              mass_star;

      for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
      {
         cell_tot_evap -= parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt * parcel_cloud.num_drop[i_pc];
      }

      parcel_cloud.temp[i_pc] = fmin(parcel_cloud.temp_starm1[i_pc], max_critical_temperature);

      enthalpy_liquid  = 0.0;
      enthalpy_liquid1 = 0.0;
      hvap             = 0.0;
      for(int isp = 0; isp < num_parcel_species; isp++)
      {
         enthalpy_liquid += parcel_cloud.mfrac[i_pc * num_parcel_species + isp] *
                            CONVERGE_table_lookup(h_table[isp], parcel_cloud.temp_tm1[i_pc]);
         enthalpy_liquid1 += parcel_cloud.mfrac[i_pc * num_parcel_species + isp] *
                             CONVERGE_table_lookup(h_table[isp], parcel_cloud.temp[i_pc]);
         // Urea model to be implemented
         hvap -= CONVERGE_table_lookup(hvap_table[isp], parcel_cloud.temp[i_pc]) *
                 parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt * parcel_cloud.num_drop[i_pc];
      }

      enth = 0.0;
      local_sie = 0.0;
      if(spray_evap_flag > 0)
      {
         if(drop_evap_source_flag == 0)
         {
            enth = CONVERGE_table_lookup(evap_species_sensible_h_table[local_evap_species_index], parcel_cloud.temp_tm1[i_pc]) * tot_evap_mass;
            local_sie = (CONVERGE_table_lookup(evap_species_h_table[local_evap_species_index], parcel_cloud.temp_tm1[i_pc]) -
                          gas_constant / CONVERGE_species_mw(sp, local_evap_species_index) * parcel_cloud.temp_tm1[i_pc]) *
                         tot_evap_mass;
         }
         else
         {
            for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
            {
               // Urea model to be implemented
               CONVERGE_index_t isp1 = CONVERGE_parcel_evap_species_lookup(isp);
               enth += CONVERGE_table_lookup(evap_species_sensible_h_table[isp1], parcel_cloud.temp_tm1[i_pc]) *
                       (-parcel_cloud.dm_dt[i_pc*num_parcel_species+isp]) * dt * parcel_cloud.num_drop[i_pc];
               local_sie += CONVERGE_table_lookup(evap_species_h_table[isp1], parcel_cloud.temp_tm1[i_pc]) -
                       gas_constant / CONVERGE_species_mw(sp, isp1) * parcel_cloud.temp_tm1[i_pc] *
                       (-parcel_cloud.dm_dt[i_pc*num_parcel_species+isp] * dt * parcel_cloud.num_drop[i_pc]);
            }
         }
      }

      //  heat of vaporization should be added to the evaporated liquid
      denth += (enth - hvap) - mass_star * (enthalpy_liquid1 - enthalpy_liquid);
      if(CONVERGE_get_int("hidden.energy_output_flag"))
      {
         region_energy_parcel_phase_change[global_region_index[node_index]] += hvap;
         region_energy_parcel_heating[global_region_index[node_index]] += mass_star*(enthalpy_liquid1-enthalpy_liquid);
         region_energy_parcel_to_gas_phase[global_region_index[node_index]] += local_sie;
      }

      /*
       * Only useful if evap_source_flag is 1 or 2.
       * Should be done for both ODE / big_evap model.
       */
      for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
      {
         evap_cell_tot_evap_species[isp] -=
            parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt * parcel_cloud.num_drop[i_pc];
         cell_tot_temp_species[isp] -= parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt *
                 parcel_cloud.num_drop[i_pc] * parcel_cloud.temp[i_pc];
      }

      // add src term for momentum
      CONVERGE_precision_t *cloud_src_evap_mom = CONVERGE_cloud_get_src_evap_mom(cloud);
      for(int ii = 0; ii < 3; ii++)
      {
         cloud_src_evap_mom[ii] = cloud_src_evap_mom[ii] + dmass * parcel_cloud.uu[i_pc][ii] / (dt * local_volume);
      }

      // update parcel radius
      radius_new[i_pc] = (radius_new[i_pc] < 1.0e-18) ? (1.0e-18) : (radius_new[i_pc]);

      parcel_cloud.radius[i_pc] = radius_new[i_pc];
   }

   *CONVERGE_cloud_get_cell_tot_evap(cloud) = cell_tot_evap;

   // adjust gas sie to account for evap.
   global_src_sie_ex[node_index] += denth / (dt * local_volume);
   global_src_unburned_enth_evap_ex[node_index] += denth / (dt * local_volume);

   // clip for too big negative src_sie_ex
   if(global_src_sie_ex[node_index] < 0.0)
   {
      // subtract min_temp so source will not result in temp below min_temp
      approx_sie   = (global_temperature[node_index] - (local_min_temp - 20.0)) * global_csubv[node_index];
      src_sie_evap = global_src_sie_ex[node_index] * dt / (global_density[node_index] + cell_tot_evap / local_volume);
      if(src_sie_evap / approx_sie < -1.00)
      {
         global_src_sie_ex[node_index]                = -1.00 * approx_sie * global_density[node_index] / dt;
         global_src_unburned_enth_evap_ex[node_index] = -1.00 * approx_sie * global_density[node_index] / dt;
      }
   }

   // adjust species_density to account for evap.

   global_src_ex_density[node_index] += cell_tot_evap / (dt * local_volume);

   if(spray_evap_flag > 0)
   {
      if(drop_evap_source_flag == 0)
      {
         local_source = cell_tot_evap / (dt * local_volume);
         if(urea_flag == 0)
         {
            if(ecfm3z_flag != 1)
            {
               global_src_species_ex[node_index * num_total_species + local_evap_species_index] += local_source;
            }
            global_src_rif_fuel_evap_ex[node_index] += local_source;
            global_ecfm_temp_evap_ex[node_index] += cell_tot_evap/(dt*global_volume[node_index]);
         }
      }
      else // drop_evap_source_flag==1 or 2
      {
         for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
         {
            // Urea model to be implemented
            CONVERGE_index_t isp1 = CONVERGE_parcel_evap_species_lookup(isp);
            local_source = evap_cell_tot_evap_species[isp] / (dt * local_volume);
            if(ecfm3z_flag > 0)
            {
               global_src_rif_fuel_evap_ex[node_index] += local_source;
               global_ecfm_temp_evap_ex[node_index] += cell_tot_temp_species[isp] / (dt * global_volume[node_index]);
            }
            // ecfm3z flag == 1 to be implemented
            global_src_species_ex[node_index * num_total_species + isp1] +=
                    evap_cell_tot_evap_species[isp] / (dt * local_volume);
            // sage_pdf_flag == 1 to be implemented
         }
      }
   }

   local_density = global_density[node_index] + global_src_ex_density[node_index] * dt;

   sum = 0.0;
   for(CONVERGE_index_t isp = 0; isp < num_gas_species; isp++)
   {
      local_species[isp] =
         (global_species_massfrac_tm1[node_index * num_total_species + isp] * global_density_tm1[node_index] +
          global_src_species_ex[node_index * num_total_species + isp] * dt) /
         local_density;
      sum = sum + local_species[isp];
   }

   for(CONVERGE_index_t isp = 0; isp < num_gas_species; isp++)
   {
      local_species[isp] = local_species[isp] / sum;
   }

   local_sensible_sie = global_sensible_sie[node_index] + global_src_sie_ex[node_index] * dt / local_density;

   temp_new =
      CONVERGE_get_temp_from_sensible_sie_massfrac(local_sensible_sie, local_species, global_pressure[node_index]);

   if(temp_new < 0.9 * global_temperature[node_index] || temp_new < min_spray_temp)
   {
      temp_new = fmax(min_spray_temp, 0.9 * global_temperature[node_index]);

      local_sensible_sie =
         CONVERGE_get_sensible_sie_from_temp_massfrac(local_species, temp_new, global_pressure[node_index]);
      global_src_sie_ex[node_index] = (local_sensible_sie - global_sensible_sie[node_index]) * local_density / dt;
      CONVERGE_set_dt_recover(CONVERGE_get_dt_max() * 10.0);
   }

   // Update profiling info
   CONVERGE_precision_t end_time = CONVERGE_mpi_wtime();
   total_evap_time += end_time - start_time;
   evap_call_count++;

   // Print profiling information periodically
   if (evap_call_count % 1000 == 0) {
       CONVERGE_int_t rank;
       CONVERGE_mpi_comm_rank(&rank);
       if (rank == 0)
       {
          print_evap_profiling();
       }
   }

   free(radius_new);
   free(moles);
   free(mol_frac);
   free(local_species);
   free(vapor_mass);
   free(vapor_mass_0);
   free(parcel_mole_fraction);
   free(evap_mass_drop_0);           /* initial evaporting droplet mass */
   free(evap_mass_drop_1);           /* droplet mass after evaporation */
   free(evap_min_radius);            /* minimum evaporation dropet radius */
   free(evap_radius);                /* evaporation dropet radius */
   free(evap_cell_tot_evap_species); /* total evaporation mass for each sepcies in a cell */
   free(evap_all_flag);
   free(parcel_species_boil_flag);
   free(cell_tot_temp_species);

      // Record evaporation calculation time
      section_end = CONVERGE_mpi_wtime();
      evap_time += section_end - section_start;
   
      // Section 4: Source terms
      section_start = CONVERGE_mpi_wtime();
   if( (parcel_boil_correlation_flag==1 && spray_evap_flag!=0) || (evap_flag_flash_boiling==1) )
   {
      free(temp_boil);
   }

   if( evap_flag_flash_boiling )
   {
      free(volume_fraction);
   }

   CONVERGE_iterator_destroy(&psp_it);
   CONVERGE_iterator_destroy(&pc_it);
   if(hidden_multi_component_diffusion_flag == 1)
   {
       free(mult_ro_mass_diff);
       free(mult_sc_num);
       free(mult_sh_num);
   }
}

/**
 * @brief reset_parcel_temp_mfrac: Reset the parcel temperature and mass fractions to tm1 values
 *
 * @param parcel: parcel of interest.
 * @param temp_drop: parcel's current temp.
 * @param temp_drop_starm1: parcel's previous iteration temp.
 */
void reset_parcel_temp_mfrac(const CONVERGE_index_t parcel_idx,
                             CONVERGE_precision_t *temp_drop,
                             CONVERGE_precision_t *temp_drop_starm1)
{
   *temp_drop        = parcel_cloud.temp[parcel_idx];
   *temp_drop_starm1 = parcel_cloud.temp[parcel_idx];

   CONVERGE_precision_t *local_mfrac     = parcel_cloud.mfrac + parcel_idx * num_parcel_species;
   CONVERGE_precision_t *local_mfrac_tm1 = parcel_cloud.mfrac_tm1 + parcel_idx * num_parcel_species;
   for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
   {
      local_mfrac[isp] = local_mfrac_tm1[isp];
   }
}

void init_tables(CONVERGE_species_t species)
{
   CONVERGE_iterator_t parcel_species_it, gas_species_it;
   CONVERGE_species_parcel_iterator_create(species, &parcel_species_it);
   CONVERGE_species_gas_iterator_create(species, &gas_species_it);

   // Get the parcel species tables
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_PVAP_ID, &pvap_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_HVAP_ID, &hvap_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_CSUBP_ID, &cp_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_DENSITY_ID, &rho_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_ENTHALPY_ID, &h_table);

   // Get the gas species tables
   load_species_tables(gas_species_it, TEMPERATURE_TABLE_ENTHALPY_ID, &evap_species_h_table);
   load_species_tables(gas_species_it, TEMPERATURE_TABLE_SENSIBLE_ENTHALPY_ID, &evap_species_sensible_h_table);
   load_species_tables(gas_species_it, TEMPERATURE_TABLE_VISC_ID, &visc_table);
   load_species_tables(gas_species_it, TEMPERATURE_TABLE_COND_ID, &cond_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
   CONVERGE_iterator_destroy(&gas_species_it);
}

static void destroy_tables(CONVERGE_species_t species)
{
   CONVERGE_iterator_t parcel_species_it, gas_species_it;
   CONVERGE_species_parcel_iterator_create(species, &parcel_species_it);
   CONVERGE_species_gas_iterator_create(species, &gas_species_it);

   // Destroy local parcel tables
   unload_species_tables(parcel_species_it, &pvap_table);
   unload_species_tables(parcel_species_it, &hvap_table);
   unload_species_tables(parcel_species_it, &cp_table);
   unload_species_tables(parcel_species_it, &rho_table);
   unload_species_tables(parcel_species_it, &h_table);

   // Destroy local gas tables
   unload_species_tables(gas_species_it, &evap_species_h_table);
   unload_species_tables(gas_species_it, &evap_species_sensible_h_table);
   unload_species_tables(gas_species_it, &visc_table);
   unload_species_tables(gas_species_it, &cond_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
   CONVERGE_iterator_destroy(&gas_species_it);
}

/**
 * @brief compute_parcel_mole_fractions: convertes parcel species mass
 *                                       fractions to mole fractions
 *
 * @param variable
 * @param mass_fraction: array contaning species mass fractions
 * @param mole_fraction: array contaning species mole fractions
 *
 * @return
 */
static CONVERGE_index_t compute_parcel_mole_fractions(CONVERGE_precision_t *mass_fraction,
                                                      CONVERGE_precision_t *mole_fraction,
                                                      CONVERGE_index_t num_parcel_species,
                                                      CONVERGE_species_t sp,
                                                      CONVERGE_index_t spray_evap_flag)
{
   CONVERGE_precision_t summation_denom, molecular_weight_mix;
   CONVERGE_index_t ii;

   if(spray_evap_flag == 0)
   {
      return 0;
   }

   summation_denom      = 0.0;
   molecular_weight_mix = 0.0;
   for(ii = 0; ii < num_parcel_species; ii++)
   {
      summation_denom +=
         mass_fraction[ii] / (CONVERGE_get_parcel_species_mw(ii) + 1.0e-30);
   }
   // Average molecular weight of the multi-component parcel.
   molecular_weight_mix = 1.0 / (summation_denom + 1.0e-30);
   for(ii = 0; ii < num_parcel_species; ii++)
   {
      mole_fraction[ii] = mass_fraction[ii] * molecular_weight_mix /
                          (CONVERGE_get_parcel_species_mw(ii) + 1.0e-30);
      if(mass_fraction[ii] <= 0.0 && mole_fraction[ii] > 0.0)
      {
         CONVERGE_logger_fatal_with_mode(ANYRANK,
                                         "Species mass fraction in a parcel is 0.0 but mole fraction is %.16f.",
                                         mole_fraction[ii]);
      }
   }

   return 0;
}
