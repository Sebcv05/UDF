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
#include <globals.h>

CONVERGE_precision_t substep_tau_coeff = 0.005;
int substep_max_n = 100;




static void spray_evap_cell(CONVERGE_cloud_t cloud);
static void reset_parcel_temp_mfrac(const CONVERGE_index_t parcel_idx,
                                    CONVERGE_precision_t *temp_drop,
                                    CONVERGE_precision_t *temp_drop_starm1);
static CONVERGE_index_t compute_parcel_mole_fractions(CONVERGE_precision_t *mass_fraction,
                                                      CONVERGE_precision_t *mole_fraction,
                                                      CONVERGE_index_t num_parcel_species,
                                                      CONVERGE_species_t sp,
                                                      CONVERGE_index_t spray_evap_flag);
static CONVERGE_precision_t calculate_lk_y1_star(CONVERGE_precision_t T_drop,
                                                  CONVERGE_precision_t P_gas,
                                                  CONVERGE_precision_t P_sat,
                                                  CONVERGE_precision_t radius,
                                                  CONVERGE_precision_t drdt_prev,
                                                  CONVERGE_precision_t mu_gas,
                                                  CONVERGE_precision_t rho_liquid,
                                                  CONVERGE_precision_t Pr_gas,
                                                  CONVERGE_precision_t Sc_gas,
                                                  CONVERGE_precision_t Y_inf,
                                                  CONVERGE_precision_t M_species,
                                                  CONVERGE_precision_t M_gas_avg,
                                                  CONVERGE_index_t lk_diagnostic_flag);

// Context for implicit LK bisection solver
typedef struct {
   CONVERGE_precision_t temp1;
   CONVERGE_precision_t global_pressure;
   CONVERGE_precision_t vapor_pres;
   CONVERGE_precision_t radius;
   CONVERGE_precision_t mol_visc;
   CONVERGE_precision_t rho_liquid;
   CONVERGE_precision_t pr_num;
   CONVERGE_precision_t sc_num;
   CONVERGE_precision_t y1_inf;
   CONVERGE_precision_t m_species;
   CONVERGE_precision_t w_0;
   CONVERGE_precision_t mass_trans_coeff_geom; // mass_trans_coeff * (2*r*rho) -> independent of bsub_d
   CONVERGE_precision_t latent_heat;
   CONVERGE_precision_t evap_mass_drop_0;
   CONVERGE_precision_t dt;
   CONVERGE_precision_t mass_drop_tm1;
   CONVERGE_index_t lk_diagnostic_flag;
} LKCtx;

// Residual function for implicit solve
static CONVERGE_precision_t lk_residual_beta(CONVERGE_precision_t beta_guess, void* ctx_ptr);
// // Fuller diffusion coefficient for NH3/N2
// static CONVERGE_precision_t fuller_diffusion_coef_nh3n2(CONVERGE_precision_t T, CONVERGE_precision_t P);
// Bisection root finder
static CONVERGE_index_t solve_lk_bisection(CONVERGE_precision_t* beta_out,
                                           CONVERGE_precision_t beta_lo,
                                           CONVERGE_precision_t beta_hi,
                                           CONVERGE_precision_t tol_x,
                                           CONVERGE_precision_t tol_f,
                                           int max_iter,
                                           void* ctx);

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

static const CONVERGE_precision_t MAX_PARCEL_TEMP_DELTA = 123.0;  // 323K - 123K = 200K = min_spray_temp (thermo.dat only goes down to 200K)
static long last_temp_clamp_warn_cycle = -1;
static int lk_diagnostic_header_written = 0;

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

   // Check and display LK flags at cycle 100 to confirm they are set correctly
   if(CONVERGE_ncyc() == 100)
   {
      printf("\n\n========================================\n");
      printf("LK FLAGS CHECK AT CYCLE 100:\n");
      printf("  lk_correction_flag = %d\n", lk_correction_flag);
      printf("  lk_diagnostic_flag = %d\n", lk_diagnostic_flag);
      printf("  lk_chi_neq_min = %.4f\n", lk_chi_neq_min);
      printf("  lk_chi_neq_max = %.4f\n", lk_chi_neq_max);
      printf("========================================\n\n");
      fflush(stdout);
   }

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
   // CONVERGE_precision_t dt = CONVERGE_simulation_dt();

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
   // DEBUG: Write to file to confirm this function is called
   static int test_write_done = 0;
   if(test_write_done == 0)
   {
      FILE *fp_test = fopen("spray_evap_test.txt", "w");
      if(fp_test != NULL)
      {
         fprintf(fp_test, "spray_evap_cell() was called at ncyc=%ld\n", CONVERGE_ncyc());
         fprintf(fp_test, "lk_correction_flag = %d\n", lk_correction_flag);
         fprintf(fp_test, "lk_diagnostic_flag = %d\n", lk_diagnostic_flag);
         fclose(fp_test);
         test_write_done = 1;
      }
   }
   
   // Start timing the entire function
   
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
   CONVERGE_precision_t *latent_heat                = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   CONVERGE_precision_t *drdt_base                  = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);


   
   


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

   // Langmuir-Knudsen parameters are loaded from user_inputs.in via read_input.c (global variables)

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
   

   // This section timing will be updated at the end of the section
   
   /////////////////////////////////////////////////////////////////////

   CONVERGE_precision_t *temp_boil = NULL;
   if( (parcel_boil_correlation_flag==1 && spray_evap_flag!=0) || (evap_flag_flash_boiling==1) )
   {
      temp_boil = SAFE_calloc(num_parcel_species, CONVERGE_precision_t);
   }

   const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);


   //These must not be below the minimum temperature of NH3 in therm.dat (currently 200 K)
   const CONVERGE_precision_t min_spray_temp          = 201.0;
   const CONVERGE_precision_t min_spray_recovery_temp = 204.0;

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

   // UDF-level check for non-physical gas temperature
   if (isnan(temp_gas) || isinf(temp_gas) || temp_gas < 100.0) {
       CONVERGE_logger_err(
           "spray_evap.c: Invalid gas temperature (%.2f) in cell %ld at ncyc %ld. Skipping evaporation for all parcels in this cell.",
           temp_gas, node_index, CONVERGE_ncyc()
       );
       return; // Exit the function, skipping evap for this entire cloud/cell
   }

   // calculate prandtl number followed by reynolds number, nusselt number and sherwood number
   CONVERGE_precision_t pr_num =
      global_mol_viscosity[node_index] * global_csubp[node_index] / global_mol_cond[node_index];
      CONVERGE_precision_t user_lifetime = 0.0;
   for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(pc_it); i_pc != -1; i_pc = CONVERGE_iterator_next(pc_it))
   {
      // UDF-level check to prevent crash from non-physical parcel temperatures
      if (isnan(parcel_cloud.temp[i_pc]) || isinf(parcel_cloud.temp[i_pc]) || parcel_cloud.temp[i_pc] < 100.0) {
          CONVERGE_logger_warn(
              "spray_evap.c: Parcel %ld (cloud %ld) has invalid temperature (%.2f) at ncyc %ld before main evap loop. Skipping parcel's contribution to evap this step.",
              i_pc, parcel_cloud.cloud_index[i_pc], parcel_cloud.temp[i_pc], CONVERGE_ncyc()
          );
          continue; // Skip to the next parcel
      }
      
      // Initialize temp_drop_0 if not already set (for parcels created before this field existed)
      if(parcel_cloud.temp_drop_0[i_pc] < 1.0) {
         parcel_cloud.temp_drop_0[i_pc] = parcel_cloud.temp[i_pc];
      }
      // if((parcel_cloud.is_child[i_pc]==1 && parcel_cloud.lifetime[i_pc]<1.0e-5) || parcel_cloud.tbt[i_pc]){
      //    user_child_flag = 1;
      //    // continue;
      // }
      //printf("\n spray_evap_cell: L501, i_pc = %ld\n  ", i_pc);
      // see Borman and Ragland 1998 edition, p. 596
      tg       = (2.0 * parcel_cloud.temp[i_pc] + temp_gas) / 3.0;
      mol_visc = 0.0;
      for(int isp = 0; isp < num_gas_species; isp++)
      {
         if (visc_table == NULL) {
             CONVERGE_logger_fatal("UDF Error in spray_evap.c: The entire 'visc_table' is NULL at ncyc %ld. Table loading failed.", CONVERGE_ncyc());
             return;
         }
         if (visc_table[isp] == NULL) {
             CONVERGE_logger_err("UDF Error in spray_evap.c: 'visc_table' for species index %d is NULL at ncyc %ld. Skipping viscosity calculation for this species.", isp, CONVERGE_ncyc());
             continue;
         }
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
      //Cap max rel_velocity for heat transfer correlations
      if(parcel_cloud.rel_vel_mag[i_pc] > 1.0e2)
      {
         parcel_cloud.rel_vel_mag[i_pc] = 1.0e2;
      }

      parcel_cloud.rey_num[i_pc] = 2.0 * parcel_cloud.radius[i_pc] * global_density[node_index] *
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
         // ******************************************************************************************************//
         // Print Droplet Data to File (Text Mode with Sampling)
         // Sampling: 5% of cells, every 5 cycles to keep file size manageable
         // Uses a random parcel index within [0, num_parcels-1] to avoid always sampling the first parcel.
         if(CONVERGE_ncyc() % 5 == 0 && num_parcels > 0)
         {
             CONVERGE_precision_t user_rand = CONVERGE_random_precision();
             if(user_rand < 0.05)
             {
                 int rank;
                 CONVERGE_mpi_comm_rank(&rank);

                 // Pick a random, in-bounds parcel index for this cloud
                 CONVERGE_precision_t rnd_idx_f = CONVERGE_random_precision();
                 CONVERGE_index_t samp_idx = (CONVERGE_index_t)(rnd_idx_f * (CONVERGE_precision_t)num_parcels);
                 if(samp_idx >= num_parcels) samp_idx = num_parcels - 1;  // safety clamp
                 if(samp_idx < 0)           samp_idx = 0;

                 CONVERGE_precision_t vmag = CONVERGE_sqrt(
                     CONVERGE_square(parcel_cloud.uu[samp_idx][0]) +
                     CONVERGE_square(parcel_cloud.uu[samp_idx][1]) +
                     CONVERGE_square(parcel_cloud.uu[samp_idx][2]));

                 CONVERGE_precision_t sim_time = CONVERGE_simulation_time_sec();

                 char filename[64];
                 sprintf(filename, "Temp_Tracker_rank_%d.txt", rank);

                 // Use append mode "a" to add to the rank-specific file without overwriting
                 FILE *fp1 = fopen(filename, "a");
                 if (fp1 != NULL)
                 {
                     // Format: cloud_idx parcel_idx temp radius num_drop lifetime vmag is_child sim_time inj_time
                     fprintf(fp1, "%ld %ld %.6e %.6e %.6e %.6e %.6e %d %.6e %.6e\n",
                             parcel_cloud.cloud_index[samp_idx],
                             parcel_cloud.parcel_index[samp_idx],
                             parcel_cloud.temp[samp_idx],
                             parcel_cloud.radius[samp_idx],
                             parcel_cloud.num_drop[samp_idx],
                             parcel_cloud.lifetime[samp_idx],
                             vmag,
                             (parcel_cloud.breakup_phase[samp_idx] >= 5) ? 1 : 0,
                             sim_time,
                             parcel_cloud.time_of_injection[samp_idx]);
                     fclose(fp1);
                 }
             }
         }
         // ******************************************************************************************************//

   //******************************************************************************************************//
   //
   //       Start Implicit Solver (i.e., calculate drop temp and mass evaporated)
   //
   //******************************************************************************************************//
   const CONVERGE_precision_t local_volume = 1.0 / (1.0 / global_volume[node_index]);   //FIXVOL
CONVERGE_precision_t user_drdt = 0.0;
CONVERGE_precision_t user_parcel_temp = 0.0;
CONVERGE_precision_t user_gas_temp = 0.0;
CONVERGE_precision_t user_radius = 0.0;
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
            // //Print global pressure and temp then calculate diffusivity and print that to check units
            // printf("\n\n\nGlobal Pressure = %e",global_pressure[node_index]);
            // printf("\nTemp = %e",temp2);
            // printf("\nD0 CONVERGE = %e",spray_evap_d0_diffuse);
            // CONVERGE_precision_t converge_D =  density1 * spray_evap_d0_diffuse * pow((temp2 / 273.0), (n_diffuse_m1))/global_density[node_index];
            // printf("\nD CONVERGE = %e",converge_D);
            // CONVERGE_precision_t user_D0 = fuller_diffusion_coef_nh3n2(temp2,global_pressure[node_index]);
            // printf("\nD calc = %e\n\n\n",user_D0);




         
         int n_sub = 1;
         CONVERGE_precision_t drdt_guess = 0.0;
         for (int isp = 0; isp < num_parcel_species; isp++) {
             drdt_guess += parcel_cloud.drdt[i_pc * num_parcel_species + isp];
         }
         
         // Diagnostics storage
         CONVERGE_precision_t diag_init_drdt = drdt_guess;
         CONVERGE_precision_t diag_init_temp = parcel_cloud.temp[i_pc];
         CONVERGE_precision_t diag_final_drdt = 0.0;
         CONVERGE_precision_t diag_final_temp = 0.0;

         if (fabs(drdt_guess) > 1.0e-12) {
             CONVERGE_precision_t tau_R = parcel_cloud.radius_tm1[i_pc] / fabs(drdt_guess);
             n_sub = (int)ceil(dt / (substep_tau_coeff * tau_R));
         }
         
         // Add a thermal check: if the explicit predictor already expects a massive Delta T over this timestep, sub-step it too!
         CONVERGE_precision_t delta_T_guess = fabs(parcel_cloud.temp[i_pc] - parcel_cloud.temp_tm1[i_pc]);
         if (delta_T_guess > 2.0) {
             int n_sub_th = (int)ceil(delta_T_guess / 2.0);
             if (n_sub_th > n_sub) n_sub = n_sub_th;
         }

         if (n_sub < 1) n_sub = 1;
         if (n_sub > substep_max_n) {
             printf("\nWARNING: Max substeps reached or exceeded! cyc=%ld pc=%ld req_n_sub=%d (capped at %d)\n", 
                    (long)CONVERGE_ncyc(), (long)i_pc, n_sub, substep_max_n);
             n_sub = substep_max_n;
         }
         
         CONVERGE_precision_t dt_sub = dt / n_sub;

         // Sub-step tm1 local variables
         CONVERGE_precision_t sub_temp_tm1 = parcel_cloud.temp_tm1[i_pc];
         CONVERGE_precision_t sub_radius_tm1 = parcel_cloud.radius_tm1[i_pc];
         CONVERGE_precision_t sub_density_tm1 = parcel_cloud.density_tm1[i_pc];
         CONVERGE_precision_t sub_mfrac_tm1[100]; // assuming max 100 species
         for(int isp=0; isp<num_parcel_species; ++isp) {
             sub_mfrac_tm1[isp] = parcel_cloud.mfrac_tm1[i_pc * num_parcel_species + isp];
         }
         CONVERGE_precision_t sub_mass_tm1 = sub_density_tm1 * (4.0 / 3.0) * 3.14159265358979323846 * (sub_radius_tm1 * sub_radius_tm1 * sub_radius_tm1);
         
         // Accumulators for cell source sums over dt
         CONVERGE_precision_t acc_dm_dt[100];
         for(int isp=0; isp<num_parcel_species; ++isp) acc_dm_dt[isp] = 0.0;

         // Start sub-loop
         CONVERGE_precision_t orig_radius_macro = parcel_cloud.radius[i_pc];
         CONVERGE_precision_t orig_density_macro = parcel_cloud.density[i_pc];
         CONVERGE_precision_t orig_temp_macro = parcel_cloud.temp[i_pc];
         for (int sub_iter = 0; sub_iter < n_sub; ++sub_iter)
         {
         // Proxy the global array element to safely isolate sub-step shrinkage physics
         parcel_cloud.radius[i_pc]  = sub_radius_tm1;
         parcel_cloud.density[i_pc] = sub_density_tm1;
         parcel_cloud.temp[i_pc]    = sub_temp_tm1;
         
         int inner_iter_flag                 = 1;
         int inner_iter                      = 0;
         CONVERGE_precision_t inner_iter_tol = 1.0;

         CONVERGE_precision_t tdrop        = (sub_iter == 0) ? parcel_cloud.temp[i_pc] : sub_temp_tm1;
         CONVERGE_precision_t tdrop_starm1 = (sub_iter == 0) ? parcel_cloud.temp[i_pc] : sub_temp_tm1;

         CONVERGE_precision_t temp_prev_timestep = sub_temp_tm1;
         if(temp_prev_timestep < min_spray_temp || temp_prev_timestep > max_critical_temperature)
         {
            CONVERGE_logger_warn(
               "spray_evap.c: Previous timestep temperature reset for parcel %ld (cloud %ld) at ncyc %ld. "
               "temp_tm1=%.3e outside [%g, %g]",
               i_pc,
               parcel_cloud.cloud_index[i_pc],
               CONVERGE_ncyc(),
               temp_prev_timestep,
               min_spray_temp,
               max_critical_temperature);

            temp_prev_timestep            = fmin(fmax(temp_prev_timestep, min_spray_temp), max_critical_temperature);
            sub_temp_tm1   = temp_prev_timestep;
            tdrop                        = temp_prev_timestep;
            tdrop_starm1                 = temp_prev_timestep;
         }

         if(parcel_cloud.temp[i_pc] > parcel_cloud.temp_drop_0[i_pc] + 2.0 || parcel_cloud.temp_tm1[i_pc] > parcel_cloud.temp_drop_0[i_pc] + 2.0)
         {
            printf("spray_evap: temp=%.3f temp_tm1=%.3f temp_drop_0+2K=%.3f radius=%.3e num_drop=%lld breakup_phase=%d lifetime=%.3e\n",
                   parcel_cloud.temp[i_pc],
                   parcel_cloud.temp_tm1[i_pc],
                   parcel_cloud.temp_drop_0[i_pc] + 2.0,
                   parcel_cloud.radius[i_pc],
                   (long long)parcel_cloud.num_drop[i_pc],
                   (int)parcel_cloud.breakup_phase[i_pc],
                   parcel_cloud.lifetime[i_pc]);
         }

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

            // UDF-level check to prevent crash from non-physical temperatures during solver recovery
            if (isnan(temp1) || isinf(temp1) || temp1 < 100.0 || isnan(temp2) || isinf(temp2) || temp2 < 100.0) {
                CONVERGE_logger_err(
                    "spray_evap.c: Invalid temperature in solver for parcel %ld (cloud %ld) at ncyc %ld. T_drop=%.2f, T_gas=%.2f. Breaking inner loop.",
                    i_pc, parcel_cloud.cloud_index[i_pc], CONVERGE_ncyc(), tdrop, temp_gas
                );
                // Reset parcel temp to a safe value from previous timestep and exit this parcel's evap calculation for this step
                parcel_cloud.temp[i_pc] = parcel_cloud.temp_tm1[i_pc];
                tdrop = parcel_cloud.temp_tm1[i_pc];
                break; // Exit the 'while' loop for this parcel
            }

            CONVERGE_precision_t thermal_cond_gas = 0.0;
            mol_visc                              = 0.0;
            for(int isp = 0; isp < num_gas_species; isp++)
            {
               if (cond_table == NULL || cond_table[isp] == NULL) {
                   CONVERGE_logger_fatal("FATAL ERROR in spray_evap.c: cond_table for species index %d is NULL at ncyc %ld. Table loading likely failed.", isp, CONVERGE_ncyc());
                   CONVERGE_mpi_abort();
               }
               thermal_cond_gas += global_species_massfrac[node_index * num_total_species + isp] *
                                   CONVERGE_table_lookup(cond_table[isp], temp2);

               if (visc_table == NULL || visc_table[isp] == NULL) {
                   CONVERGE_logger_fatal("FATAL ERROR in spray_evap.c: visc_table for species index %d is NULL at ncyc %ld. Table loading likely failed.", isp, CONVERGE_ncyc());
                   CONVERGE_mpi_abort();
               }
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

            CONVERGE_precision_t mass_drop_tm1 = sub_mass_tm1;

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
            
            // Skip entire parcel if radius is invalid (completely evaporated or uninitialized)
            if(isnan(parcel_cloud.radius[i_pc]) || parcel_cloud.radius[i_pc] < 1.0e-10)
            {
               printf("[EVAP_SKIP] Parcel i_pc=%d has invalid radius=%.3e m at ncyc=%ld - skipping entire parcel\n",
                      (int)i_pc, parcel_cloud.radius[i_pc], CONVERGE_ncyc());
               // Set all species evaporation terms to zero
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = 0.0;
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
               }
               continue;  // Skip to next parcel
            }
            
            for(CONVERGE_index_t isp = 0, isp1 = CONVERGE_iterator_first(psp_it); isp < num_parcel_species;
                isp++, isp1                    = CONVERGE_iterator_next(psp_it))
            {
               evap_all_flag[isp] = 0;
               if(sub_mfrac_tm1[isp] <= 0.0)
            {
               parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = 0.0;
               evap_mass_drop_0[isp]                               = 0.0;
               latent_heat[isp]                                    = 0.0;
               drdt_base[isp]                                      = parcel_cloud.drdt[i_pc * num_parcel_species + isp];
               continue;
            }

            evap_mass_drop_0[isp] = mass_drop_tm1 * sub_mfrac_tm1[isp];

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
               latent_heat[isp]                                    = 0.0;
               drdt_base[isp]                                      = parcel_cloud.drdt[i_pc * num_parcel_species + isp];
               continue;
            }

               CONVERGE_precision_t y1_star = 0.0;
               if(spray_evap_flag > 0)
               {
                  // Apply Langmuir-Knudsen correction if enabled
                  if(lk_correction_flag == 1)
                  {
                     // Get values for safety checks and debugging
                     CONVERGE_precision_t drdt_prev = parcel_cloud.drdt[i_pc * num_parcel_species + isp];
                     
                     // Print values before safety checks (always print for invalid parcels)
                     static int debug_print_count = 0;
                     int print_this = 0;
                     
                     if(debug_print_count < 10)
                     {
                        print_this = 1;
                        debug_print_count++;
                     }
                     
                     // Also print if we detect potential issues
                     if(isnan(temp1) || isinf(temp1) || temp1 < 100.0 || temp1 > 1000.0 ||
                        isnan(drdt_prev) || isinf(drdt_prev) || fabs(drdt_prev) > 1.0e6 ||
                        isnan(parcel_cloud.radius[i_pc]) || parcel_cloud.radius[i_pc] < 1.0e-10)
                     {
                        print_this = 1;
                     }
                     
                     if(print_this)
                     {
                        printf("[LK_VALUES] ncyc=%ld, i_pc=%d, isp=%d:\n", CONVERGE_ncyc(), (int)i_pc, (int)isp);
                        printf("  temp1=%.3f K, vapor_pres=%.1f Pa, P_gas=%.1f Pa\n", 
                               temp1, vapor_pres, global_pressure[node_index]);
                        printf("  radius=%.3e m, drdt_prev=%.3e m/s\n", 
                               parcel_cloud.radius[i_pc], drdt_prev);
                        printf("  breakup_phase=%d, lifetime=%.3e s\n", 
                               parcel_cloud.breakup_phase[i_pc], parcel_cloud.lifetime[i_pc]);
                        printf("  isnan(temp1)=%d, isinf(temp1)=%d, isnan(drdt)=%d, isinf(drdt)=%d\n",
                               isnan(temp1), isinf(temp1), isnan(drdt_prev), isinf(drdt_prev));
                     }
                     
                     // Safety checks before calling LK model
                     int skip_lk = 0;
                     
                     // Check temperature is valid
                     if(isnan(temp1) || isinf(temp1) || temp1 < 100.0 || temp1 > 1000.0)
                     {
                        printf("[LK_SKIP] Invalid temp1=%.3e at ncyc=%ld, i_pc=%d\n", 
                               temp1, CONVERGE_ncyc(), (int)i_pc);
                        skip_lk = 1;
                     }
                     
                     // Check drdt is valid (not NaN, Inf, or unreasonably large)
                     if(isnan(drdt_prev) || isinf(drdt_prev) || fabs(drdt_prev) > 1.0e6)
                     {
                        printf("[LK_SKIP] Invalid drdt_prev=%.3e at ncyc=%ld, i_pc=%d\n", 
                               drdt_prev, CONVERGE_ncyc(), (int)i_pc);
                        skip_lk = 1;
                     }
                     
                     // Skip LK for child parcels in their first timesteps (drdt not yet established)
                     if(parcel_cloud.breakup_phase[i_pc] >= 5 && parcel_cloud.lifetime[i_pc] < 1.0e-6)
                     {
                        if(debug_print_count < 10)
                        {
                           printf("[LK_SKIP] Child parcel i_pc=%d, lifetime=%.3e s < 1e-6 at ncyc=%ld\n", 
                                  (int)i_pc, parcel_cloud.lifetime[i_pc], CONVERGE_ncyc());
                        }
                        skip_lk = 1;
                     }
                     
                     if(skip_lk == 1)
                     {
                        // Fall back to classical model
                        if(drop_evap_source_flag == 0)
                        {
                           y1_star = CONVERGE_species_mw(sp, local_evap_species_index) /
                                     (CONVERGE_species_mw(sp, local_evap_species_index) +
                                      (w_0 * (global_pressure[node_index] / vapor_pres - 1.0)));
                        }
                        else
                        {
                           y1_star = parcel_mole_fraction[isp] * vapor_pres /
                                   global_pressure[node_index] * CONVERGE_get_parcel_species_mw(isp) / w_0_denom;
                        }
                     }
                     else
                     {
                        // Get liquid density for this species with safety check
                        if(rho_table == NULL || rho_table[isp] == NULL)
                        {
                           printf("[LK_ERROR] rho_table not initialized for species %d\n", (int)isp);
                           y1_star = 1.0e-10;
                           continue;
                        }
                        
                        CONVERGE_precision_t rho_liquid = CONVERGE_table_lookup(rho_table[isp], temp1);

                        // NEW: Implicit LK solution for small droplets (D < 50 um)
                        const CONVERGE_precision_t D_threshold = 50.0e-6; // 50 microns
                        CONVERGE_precision_t D_drop = parcel_cloud.radius[i_pc] * 2.0;

                        // Pre-calculate Sherwood / mass trans coeff terms for both approaches
                        CONVERGE_precision_t mass_trans_coeff_geom = 0.0;
                        if(hidden_multi_component_diffusion_flag == 1)
                        {
                            mass_trans_coeff_geom = spray_scale_mass_trans_coeff_spray * mult_sh_num[i_pc * num_parcel_species + isp] * (mol_visc / mult_sc_num[isp]);
                        }
                        else
                        {
                           mass_trans_coeff_geom = spray_scale_mass_trans_coeff_spray * parcel_cloud.v_sh[i_pc] * (mol_visc / sc_num);
                        }
                        
                        // Decide between implicit and explicit LK
                        int implicit_success = 0;
                        if (D_drop < D_threshold) 
                        {
                           // Setup context for the Residual Evaluator
                           LKCtx lk_ctx;
                           lk_ctx.temp1 = temp1;
                           lk_ctx.global_pressure = global_pressure[node_index];
                           lk_ctx.vapor_pres = vapor_pres;
                           lk_ctx.radius = parcel_cloud.radius[i_pc];
                           lk_ctx.mol_visc = mol_visc;
                           lk_ctx.rho_liquid = rho_liquid;
                           lk_ctx.pr_num = pr_num;
                           lk_ctx.sc_num = sc_num;
                           lk_ctx.y1_inf = vapor_mass[isp] / mass;
                           lk_ctx.m_species = CONVERGE_get_parcel_species_mw(isp);
                           lk_ctx.w_0 = w_0;
                           lk_ctx.mass_trans_coeff_geom = mass_trans_coeff_geom;
                           lk_ctx.latent_heat = hvap;
                           lk_ctx.evap_mass_drop_0 = evap_mass_drop_0[isp];
                           lk_ctx.dt = dt_sub;
                           lk_ctx.mass_drop_tm1 = mass_drop_tm1;
                           lk_ctx.lk_diagnostic_flag = 0; // Turn off inner loop diagnostics
                           
                           // Guess Beta from explicit evaluation as the upper bound
                           CONVERGE_precision_t C_beta_guess = (rho_liquid * pr_num * parcel_cloud.radius[i_pc]) / (mol_visc + 1.0e-30);
                           CONVERGE_precision_t beta_guess_hi = -C_beta_guess * drdt_prev;
                           if (beta_guess_hi < 1.0e-6) { beta_guess_hi = 1.0e-6; }
                           
                           CONVERGE_precision_t beta_lo = 0.0;
                           CONVERGE_precision_t beta_root = 0.0;
                           
                           CONVERGE_index_t slv_status = solve_lk_bisection(&beta_root, beta_lo, beta_guess_hi, 1e-8, 1e-10, 40, &lk_ctx);
                           
                           if (slv_status == 0) // successfully converged
                           {
                              // Transform back to drdt to lock final y1_star
                              CONVERGE_precision_t drdt_root = -beta_root / (C_beta_guess + 1.0e-30);
                              
                              y1_star = calculate_lk_y1_star(
                                 temp1, global_pressure[node_index], vapor_pres, parcel_cloud.radius[i_pc],
                                 drdt_root, mol_visc, rho_liquid, pr_num, sc_num, vapor_mass[isp] / mass,
                                 CONVERGE_get_parcel_species_mw(isp), w_0, lk_diagnostic_flag
                              );

                              implicit_success = 1;
                              
                              static int implicit_dbg_count = 0;
                              if (implicit_dbg_count < 10) {
                                  int rank;
                                  CONVERGE_mpi_comm_rank(&rank);
                                  printf("[Rank %d LK_IMPLICIT] Success. i_pc=%d isp=%d. beta_guess=%.3e beta_root=%.3e drdt_prev=%.3e drdt_root=%.3e Ys=%.6f\n",
                                         rank, (int)i_pc, (int)isp, beta_guess_hi, beta_root, drdt_prev, drdt_root, y1_star);
                                  implicit_dbg_count++;
                              }
                           }
                           else 
                           {
                              static int fail_dbg_count = 0;
                              if (fail_dbg_count < 20) {
                                  int rank;
                                  CONVERGE_mpi_comm_rank(&rank);
                                  printf("[Rank %d LK_IMPLICIT_FAIL] Code %d. i_pc=%d isp=%d. Falling back to explicit LK.\n", 
                                         rank, (int)slv_status, (int)i_pc, (int)isp);
                                  fail_dbg_count++;
                              }
                           }
                        }

                        if (!implicit_success)
                        {
                           // Use classical explicit LK model if we failed, or if diameter >= 50 microns
                           y1_star = calculate_lk_y1_star(
                              temp1,                                           // T_drop
                              global_pressure[node_index],                     // P_gas
                              vapor_pres,                                      // P_sat
                              parcel_cloud.radius[i_pc],                       // radius
                              drdt_prev,                                       // drdt_prev
                              mol_visc,                                        // mu_gas
                              rho_liquid,                                      // rho_liquid
                              pr_num,                                          // Pr_gas
                              sc_num,                                          // Sc_gas
                              vapor_mass[isp] / mass,                          // Y_inf
                              CONVERGE_get_parcel_species_mw(isp),            // M_species
                              w_0,                                             // M_gas_avg
                              lk_diagnostic_flag                               // diagnostic flag
                           );
                        }
                     }
                  }
                  else
                  {
                     // Classical model (existing calculation)
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


               //USER MAXWELL MODEL IMPLEMENTATION
               //Sherwood number (modidfied)
               //   CONVERGE_precision_t user_sh = 2.0 + 0.6 * sqrt(parcel_cloud.rey_num[i_pc]) * (CONVERGE_cbrt(sc_num));

               // mass_trans_coeff = spray_scale_mass_trans_coeff_spray * user_sh * (mol_visc / sc_num) /
                                

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

              //USER MAXWELL MODEL IMPLEMENTATION
//        if (evap_flag_flash_boiling == 1)
// {
//     double super_heat_degree = tdrop - temp_boil[isp];
//     if (super_heat_degree > 0.2)
//     {
//         // --- molar masses (kg/mol) ---
//         const CONVERGE_precision_t M_NH3 = 17.031e-3;
//         const CONVERGE_precision_t M_N2  = 28.014e-3;

//         // --- diffusivity (verify mol_visc is kinematic viscosity [m^2/s]) ---
//         CONVERGE_precision_t user_D = mol_visc / sc_num;
//          // printf("\n D = %e , T= %f K, Sh = %e\n", user_D, temp1, parcel_cloud.v_sh[i_pc]);
//         // --- universal gas constant (J/mol/K) ---
//         const CONVERGE_precision_t R_univ = 8.314462618;

//         // --- convert mass fraction y1 (NH3 mass fraction) -> mole fraction x_NH3 ---
//         CONVERGE_precision_t denom_x1 = (M_N2 * y1) + (M_NH3 * (1.0 - y1));
//         if (denom_x1 <= 0.0) denom_x1 = 1e-30;
//         CONVERGE_precision_t user_x1 = (M_N2 * y1) / denom_x1;
//         if (user_x1 < 0.0) user_x1 = 0.0;
//         else if (user_x1 > 1.0) user_x1 = 1.0;

//         // --- partial pressure of NH3 in the cell [Pa] ---
//         CONVERGE_precision_t P_tot = global_pressure[node_index]; // MUST be Pa
//         CONVERGE_precision_t user_Pv = user_x1 * P_tot;

//         // --- saturation vapor pressure of ammonia at droplet surface temperature [Pa] ---
//         CONVERGE_precision_t P_sat = CONVERGE_table_lookup(pvap_table[isp], temp1); // Pa

    

//         // --- compute x_star if you still need it (kept for reference) ---
//         CONVERGE_precision_t denom_xstar = (M_N2 * y1_star) + (M_NH3 * (1.0 - y1_star));
//         if (denom_xstar <= 0.0) denom_xstar = 1e-30;
//         CONVERGE_precision_t user_xstar = (M_N2 * y1_star) / denom_xstar;
//         if (user_xstar < 0.0) user_xstar = 0.0;
//         else if (user_xstar > 1.0) user_xstar = 1.0;

//          //Psat (LK Correction with user_xstar)
//          // CONVERGE_precision_t user_Ps = user_xstar * P_sat;
//          CONVERGE_precision_t user_Ps = P_sat;
//         // --- denominator for dr/dt: r * R * rho_liquid  (rho must be liquid density [kg/m^3]) ---
//         CONVERGE_precision_t r = parcel_cloud.radius[i_pc];
//         if (r < 1e-12) r = 1e-12;
//         CONVERGE_precision_t rho_l = parcel_cloud.density[i_pc];
//         if (rho_l <= 0.0) rho_l = 1.0; // better to assert/throw in real code

//         CONVERGE_precision_t user_denom = r * R_univ * rho_l;

//         // --- concentration difference term (mol/m^3) ---
//         CONVERGE_precision_t Td = temp1;      // droplet surface temp [K]
//         CONVERGE_precision_t Ta = temp_gas;   // gas temp [K]
//         if (Td <= 0.0) Td = 1e-6;
//         if (Ta <= 0.0) Ta = 1e-6;

//         CONVERGE_precision_t conc_term = (user_Ps / (R_univ * Td)) - (user_Pv / (R_univ * Ta));

//         // --- final dr/dt (m/s) ---
//         CONVERGE_precision_t Sh = parcel_cloud.v_sh[i_pc]; // must be dimensionless
//         CONVERGE_precision_t drdt = - Sh * user_D * M_NH3 * conc_term / user_denom;

//         parcel_cloud.drdt[i_pc * num_parcel_species + isp] = drdt;
//     }
// }

               //printf("\n spray_evap_cell L815 ");
               //Zero for first 1e-6 s of child's lifetime to improve stability 
               if(parcel_cloud.breakup_phase[i_pc] >= 5 && parcel_cloud.lifetime[i_pc] < 1.0e-6)
               {
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] = -1.0e-7;
               }


               if(parcel_cloud.is_child[i_pc]==2) // not set to 2 so inactive
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
            //Monitor radius chnage rate
            user_drdt = parcel_cloud.drdt[i_pc * num_parcel_species + isp];
            user_parcel_temp = parcel_cloud.temp[i_pc];
            user_gas_temp = temp_gas;
            user_radius = parcel_cloud.radius[i_pc];
            user_lifetime = parcel_cloud.lifetime[i_pc];
            if(parcel_cloud.radius[i_pc] < evap_min_radius[isp])
            {
               parcel_cloud.drdt[i_pc * num_parcel_species + isp] =
                  -((parcel_cloud.radius[i_pc] - evap_min_radius[isp]) / dt_sub);
                  evap_all_flag[isp] = 1;
               }

      //                   // Print  Droplet Data to File 5% sample rate to see drdt values 
      //                   CONVERGE_precision_t user_rand1 = CONVERGE_random_precision();
      //                   if(user_rand1<0.05){
      //                   char *filename1 = "Evap_Tracker.txt";
      //                   FILE *fp1 = fopen("Evap_Tracker.txt", "a");
      //                   if (fp1 == NULL)
      //                   {
      //                      printf("Error opening the file %s", filename1);
      //                   }
      //                   fprintf(fp1, "%i    %i    %f    %e    %e    %e\n", parcel_cloud.cloud_index[0], parcel_cloud.parcel_index[0], parcel_cloud.temp[0], parcel_cloud.radius[0], parcel_cloud.lifetime[0],parcel_cloud.drdt[i_pc * num_parcel_species + isp]);
      //                   fclose(fp1);
      //                   // ******************************************************************************************************//
      // }



               //Cap maximum rate of radius change  - user routine
            

               // if(parcel_cloud.drdt[i_pc * num_parcel_species + isp] <-1.0e-2){
               //    parcel_cloud.drdt[i_pc * num_parcel_species + isp] = -1.0e-2;
               // }

               //  again don't allow condensation
               if(parcel_cloud.drdt[i_pc * num_parcel_species + isp] >= 0.0)
               {
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
               }
               
               // Prevent evaporation for parent parcels (pre-breakup)
               // Thermal breakup takes precedence over evaporation
               // if(parcel_cloud.is_child[i_pc] == 0)
               // {
               //    parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
               // }

               // temperature can't be above critical temp
               const CONVERGE_precision_t isp_tcrit = CONVERGE_species_tcrit(sp, isp1);
               if(tdrop >= (isp_tcrit - 1.0e-6))   // make the species disappear
               {
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] =
                     -((parcel_cloud.radius[i_pc] - evap_min_radius[isp]) / dt_sub);
                  evap_all_flag[isp] = 1;
               }

               evap_radius[isp] = parcel_cloud.radius[i_pc] + (parcel_cloud.drdt[i_pc * num_parcel_species + isp] * dt_sub);
               evap_mass_drop_1[isp] =
                  (4.0 / 3.0) * PI * parcel_cloud.density[i_pc] * (CONVERGE_cube(evap_radius[isp]));
               parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = (evap_mass_drop_1[isp] - mass_drop_tm1) / dt_sub;

               // Do not evaporate more mass than available
               if((-parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub) > evap_mass_drop_0[isp])
               {
                  parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = -evap_mass_drop_0[isp] / dt_sub;
               }
            drdt_base[isp] = parcel_cloud.drdt[i_pc * num_parcel_species + isp];
            latent_heat[isp] = hvap;
            vaporization_term =
               vaporization_term + (parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub) * hvap;
            }

            mass_drop_new = mass_drop_tm1;
            for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
            {
               mass_drop_new += parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub;
            }

            if(mass_drop_new > 1.0e-36) /* normal update */
            {
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  parcel_cloud.mfrac[i_pc * num_parcel_species + isp] =
                     (evap_mass_drop_0[isp] + parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub) /
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
                     sub_mfrac_tm1[isp];
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
                     sub_mfrac_tm1[isp];
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
                     (bsub_d_avg == 0.0) ? 0.0 : dt_sub * drop_area * heat_trans_coeff * log(1.0 + bsub_d_avg) / bsub_d_avg;
                     if(evap_flag_flash_boiling==1)
                     {
                        double superheatdegreee = tdrop - temp_gas;
                        if(superheatdegreee>0)
                        {
                                             (bsub_d_avg == 0.0) ? 0.0 : dt_sub * drop_area * heat_trans_coeff ;
                        }
                     }
               }
               else if(spray_evap_flag == 2)
               {
                  cond_term1 = dt_sub * drop_area_1 * heat_trans_coeff * (1.0 / (pow((1.0 + bsub_d_avg), 0.678)));
               }
               else if(spray_evap_flag == 0)
               {
                  cond_term1 = dt_sub * drop_area * heat_trans_coeff;
               }
            }
            //Turn of spalding number correlation for children after breakup 
           

            CONVERGE_precision_t denom = cond_term1 + (csubp_liquid * mass_drop_new);

            tdrop = (denom == 0.0) ? tdrop_starm1
                                   : ((csubp_liquid * mass_drop_new * parcel_cloud.temp[i_pc]) + vaporization_term +
                                      (cond_term1 * temp_gas) + spray_wall_heat_source) /
                                        denom;

            tdrop = omega * (tdrop - tdrop_starm1) + tdrop_starm1;

            const CONVERGE_precision_t delta_temp_candidate = tdrop - temp_prev_timestep;
            const CONVERGE_precision_t min_allowed =
               fmax(min_spray_temp, temp_prev_timestep - MAX_PARCEL_TEMP_DELTA);
            const CONVERGE_precision_t max_allowed =
               fmin(max_critical_temperature, temp_prev_timestep + MAX_PARCEL_TEMP_DELTA);
            CONVERGE_precision_t clamp_target_temp = tdrop;
            if(clamp_target_temp < min_allowed)
            {
               clamp_target_temp = min_allowed;
            }
            else if(clamp_target_temp > max_allowed)
            {
               clamp_target_temp = max_allowed;
            }

            const CONVERGE_precision_t target_delta = clamp_target_temp - temp_prev_timestep;
            if(fabs(delta_temp_candidate - target_delta) > 1.0e-12)
            {
               CONVERGE_precision_t clamp_ratio = 0.0;
               if(fabs(delta_temp_candidate) > 1.0e-14)
               {
                  clamp_ratio = target_delta / delta_temp_candidate;
               }

               if(clamp_ratio < 0.0)
               {
                  clamp_ratio = 0.0;
               }
               if(clamp_ratio > 1.0)
               {
                  clamp_ratio = 1.0;
               }

               // mark recovery so outer loop can relax
               do_recovery     = 1;
               inner_iter_flag = 1;
               const long current_cycle = CONVERGE_ncyc();
               if(last_temp_clamp_warn_cycle != current_cycle)
               {
                  last_temp_clamp_warn_cycle = current_cycle;
                  CONVERGE_logger_warn(
                     "spray_evap.c: Temperature clamp applied for parcel %ld (cloud %ld) at ncyc %ld. "
                     "DeltaT request=%.2f -> %.2f K, clamp_ratio=%.3f, temp = %.2f K, temp_tm1 = %.2f K, radius = %.2e, breakup_phase = %d",
                     i_pc,
                     parcel_cloud.cloud_index[i_pc],
                     current_cycle,
                     delta_temp_candidate,
                     target_delta,
                     clamp_ratio,
                     parcel_cloud.temp[i_pc],
                     parcel_cloud.temp_tm1[i_pc],
                     parcel_cloud.radius[i_pc],
                     parcel_cloud.breakup_phase[i_pc]);
               }

               CONVERGE_precision_t *local_mfrac =
                  parcel_cloud.mfrac + i_pc * num_parcel_species;

               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  local_mfrac[isp] = sub_mfrac_tm1[isp];
               }

               vaporization_term = 0.0;
               mass_drop_new     = mass_drop_tm1;

               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  parcel_cloud.drdt[i_pc * num_parcel_species + isp] =
                     drdt_base[isp] * clamp_ratio;

                  evap_all_flag[isp] = 0;

                  if(parcel_cloud.drdt[i_pc * num_parcel_species + isp] >= 0.0)
                  {
                     parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
                  }

                  if(parcel_cloud.radius[i_pc] < evap_min_radius[isp])
                  {
                     parcel_cloud.drdt[i_pc * num_parcel_species + isp] =
                        -((parcel_cloud.radius[i_pc] - evap_min_radius[isp]) / dt_sub);
                     evap_all_flag[isp] = 1;
                  }

                  evap_radius[isp] = parcel_cloud.radius[i_pc] +
                                     (parcel_cloud.drdt[i_pc * num_parcel_species + isp] * dt_sub);
                  if(evap_radius[isp] < 0.0)
                  {
                     evap_radius[isp] = 0.0;
                  }

                  evap_mass_drop_1[isp] =
                     (4.0 / 3.0) * PI * parcel_cloud.density[i_pc] * (CONVERGE_cube(evap_radius[isp]));
                  parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] =
                     (evap_mass_drop_1[isp] - mass_drop_tm1) / dt_sub;

                  if((-parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub) > evap_mass_drop_0[isp])
                  {
                     parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = -evap_mass_drop_0[isp] / dt_sub;
                  }

                  vaporization_term +=
                     (parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub) * latent_heat[isp];
               }

               mass_drop_new = mass_drop_tm1;
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  mass_drop_new += parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub;
               }

               if(mass_drop_new > 1.0e-36) /* normal update */
               {
                  for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
                  {
                     local_mfrac[isp] =
                        (evap_mass_drop_0[isp] + parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub) /
                        (mass_drop_new);
                     if(local_mfrac[isp] < 0.0 || evap_all_flag[isp] == 1)
                     {
                        local_mfrac[isp] = 0.0;
                     }
                  }
               }
               else /* if mass drop less than 1.0e-36 use old mass fractions */
               {
                  mass_drop_new = 1.0e-36;
                  for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
                  {
                     local_mfrac[isp] = sub_mfrac_tm1[isp];
                  }
               }

               sum = 0.0;
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  sum += local_mfrac[isp];
               }
               if(sum > 0.1) /* normalize */
               {
                  for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
                  {
                     local_mfrac[isp] = local_mfrac[isp] / sum;
                  }
               }
               else /* all species evaporated */
               {
                  for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
                  {
                     local_mfrac[isp] = sub_mfrac_tm1[isp];
                  }
               }

               CONVERGE_precision_t updated_csubp_liquid = 0.0;
               for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
               {
                  updated_csubp_liquid +=
                     CONVERGE_table_lookup(cp_table[isp], temp1) * local_mfrac[isp];
               }
               csubp_liquid = updated_csubp_liquid;

               radius_new[i_pc] =
                  (3.0 / 4.0) * mass_drop_new / (PI * parcel_cloud.density[i_pc]);
               radius_new[i_pc] = (radius_new[i_pc] > 0.0) ? CONVERGE_cbrt(radius_new[i_pc]) : 0.0;
               radius_new[i_pc] =
                  fmin(radius_new[i_pc], parcel_cloud.radius[i_pc]);

               const CONVERGE_precision_t denom =
                  cond_term1 + (csubp_liquid * mass_drop_new);

               CONVERGE_precision_t tdrop_local = tdrop_starm1;
               if(denom != 0.0)
               {
                  tdrop_local = ((csubp_liquid * mass_drop_new * sub_temp_tm1) + vaporization_term +
                                 (cond_term1 * temp_gas) + spray_wall_heat_source) /
                                denom;
               }

               tdrop = omega * (tdrop_local - tdrop_starm1) + tdrop_starm1;

               if(tdrop < min_allowed)
               {
                  tdrop = min_allowed;
               }
               else if(tdrop > max_allowed)
               {
                  tdrop = max_allowed;
               }
            }

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
               if(omega < (CONVERGE_precision_t)(1.0 / 16.0))
               {
                  CONVERGE_logger_warn("spray_evap.c: RECOVERY in inner loop for parcel %ld, cloud %ld at ncyc %ld. tdrop=%.2f, temp_gas=%.2f, breakup_phase=%d, radius=%e, radius_tm1=%e",
                                       i_pc,
                                       parcel_cloud.cloud_index[i_pc],
                                       CONVERGE_ncyc(),
                                       tdrop,
                                       temp_gas,
                                       parcel_cloud.breakup_phase[i_pc],
                                       parcel_cloud.radius[i_pc],
                                       parcel_cloud.radius_tm1[i_pc]);
                  CONVERGE_logger_warn("spray_evap.c: RECOVERY diagnostics: tdrop_starm1=%.4f, omega=%.4f, csubp_liquid=%.4e, mass_drop_new=%.4e, mass_drop_tm1=%.4e, vaporization_term=%.4e, cond_term1=%.4e, denom=%.4e",
                                       tdrop_starm1,
                                       omega,
                                       csubp_liquid,
                                       mass_drop_new,
                                       mass_drop_tm1,
                                       vaporization_term,
                                       cond_term1,
                                       denom);
               }
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
            if(fabs(tdrop - sub_temp_tm1) > 10)
            {
              printf("\nDELTA TEMP > 10 K (sub_iter=%d), cyc=%ld pc=%ld n_sub=%d, radius=%.3e, rey_num=%.3e, tdrop = %e, sub_temp_tm1 = %e\n", 
                     sub_iter, (long)CONVERGE_ncyc(), (long)i_pc, n_sub, parcel_cloud.radius[i_pc], parcel_cloud.rey_num[i_pc], tdrop, sub_temp_tm1);
            }
            tdrop_starm1 = tdrop;
         }
         
            // Accumulate dm_dt
            for(int isp=0; isp<num_parcel_species; isp++) {
                acc_dm_dt[isp] += parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub;
            }
            
            diag_final_drdt = 0.0;
            for(int isp=0; isp<num_parcel_species; isp++) {
                diag_final_drdt += parcel_cloud.drdt[i_pc * num_parcel_species + isp];
            }
            diag_final_temp = tdrop;

            // Update sub_tm1 variables
            sub_temp_tm1 = tdrop;
            sub_radius_tm1 = radius_new[i_pc];
            sub_density_tm1 = parcel_cloud.density[i_pc];
            sub_mass_tm1 = sub_density_tm1 * (4.0 / 3.0) * 3.14159265358979323846 * (sub_radius_tm1 * sub_radius_tm1 * sub_radius_tm1);
            for(int isp=0; isp<num_parcel_species; ++isp) {
                sub_mfrac_tm1[isp] = parcel_cloud.mfrac[i_pc * num_parcel_species + isp];
            }
         } // end of sub_iter loop
         
         // Restore macro properties to maintain global gas CFD iteration safety
         parcel_cloud.radius[i_pc]  = orig_radius_macro;
         parcel_cloud.density[i_pc] = orig_density_macro;
         parcel_cloud.temp[i_pc]    = orig_temp_macro;
         
         CONVERGE_precision_t delta_T_check = diag_final_temp - diag_init_temp;
         if (n_sub > 1 && fabs(delta_T_check) > 10.0){
             printf("SUBSTEP_ACTIVATED: cyc=%ld pc=%ld n_sub=%d init_T=%.2f final_T=%.2f DELTA = %.2f init_drdt=%.3e final_drdt=%.3e\n",
                    (long)CONVERGE_ncyc(), (long)i_pc, n_sub, diag_init_temp, diag_final_temp, delta_T_check, diag_init_drdt, diag_final_drdt);
         }

         for(int isp=0; isp<num_parcel_species; isp++) {
             parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = acc_dm_dt[isp] / dt;
         }

         parcel_cloud.temp_starm1[i_pc] = sub_temp_tm1;
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

      //Print final radius change rate
   // if(user_child_flag==1){

   // if( user_lifetime >1.0e-4)
   // {
   //    printf("\n spray_evap_cell: L1310, child? %d radius = %e, temperature = %f, gas temperature = %f, user_drdt = %e\n  ", user_child_flag, user_radius, user_parcel_temp, user_gas_temp, user_drdt);
   // }
   // }

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
      if (radius_new[i_pc] < 1.0e-18) {
         static int clamp_count = 0;
         if (clamp_count < 5) {
            printf("[EVAP_CLAMP] p_idx=%ld, radius_new=%.3e < 1e-18, clamping to 1e-18, breakup_phase=%d\n",
                   i_pc, radius_new[i_pc], parcel_cloud.breakup_phase[i_pc]);
            clamp_count++;
         }
         radius_new[i_pc] = 1.0e-18;
      }

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
      CONVERGE_logger_warn("spray_evap.c: RECOVERY at end of cell for node %ld at ncyc %ld. temp_new=%.2f, temp_gas=%.2f",
                                       node_index, CONVERGE_ncyc(), temp_new, global_temperature[node_index]);
      temp_new = fmax(min_spray_temp, 0.9 * global_temperature[node_index]);

      local_sensible_sie =
         CONVERGE_get_sensible_sie_from_temp_massfrac(local_species, temp_new, global_pressure[node_index]);
      global_src_sie_ex[node_index] = (local_sensible_sie - global_sensible_sie[node_index]) * local_density / dt;
      CONVERGE_set_dt_recover(CONVERGE_get_dt_max() * 10.0);
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
   free(latent_heat);
   free(drdt_base);

      
      
      
      // ;
      // printf("\n==========================");
      // printf("\n spray_evap.c total time = %e ms\n\n",total_time*1000);
      // printf("\ninit_time_frac = %f \%\n",init_time_frac*100);
      // printf("\nboil_time_frac = %f \%\n",boil_time_frac*100);
      // printf("\nevap_time_frac = %f \%\n",evap_time_frac*100);
      // printf("\nsource_time_frac = %f \%\n",source_time_frac*100);
      // printf("\n==========================");


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
 * @brief calculate_lk_y1_star: Calculate surface vapor mass fraction using Langmuir-Knudsen correction
 *
 * @param T_drop: Droplet temperature [K]
 * @param P_gas: Gas pressure [Pa]
 * @param P_sat: Saturation pressure at droplet temperature [Pa]
 * @param radius: Droplet radius [m]
 * @param drdt_prev: Rate of radius change from previous iteration [m/s]
 * @param mu_gas: Gas dynamic viscosity [Pa·s]
 * @param rho_liquid: Liquid density [kg/m³]
 * @param Pr_gas: Prandtl number (gas) [-]
 * @param Sc_gas: Schmidt number (gas) [-]
 * @param Y_inf: Ambient vapor mass fraction [-]
 * @param M_species: Molecular weight of evaporating species [kg/kmol]
 * @param M_gas_avg: Average molecular weight of gas [kg/kmol]
 * @param lk_diagnostic_flag: Flag to enable diagnostic printing
 *
 * @return Surface vapor mass fraction Y_v_s corrected by L-K model
 */
CONVERGE_precision_t calculate_lk_y1_star(CONVERGE_precision_t T_drop,
                                          CONVERGE_precision_t P_gas,
                                          CONVERGE_precision_t P_sat,
                                          CONVERGE_precision_t radius,
                                          CONVERGE_precision_t drdt_prev,
                                          CONVERGE_precision_t mu_gas,
                                          CONVERGE_precision_t rho_liquid,
                                          CONVERGE_precision_t Pr_gas,
                                          CONVERGE_precision_t Sc_gas,
                                          CONVERGE_precision_t Y_inf,
                                          CONVERGE_precision_t M_species,
                                          CONVERGE_precision_t M_gas_avg,
                                          CONVERGE_index_t lk_diagnostic_flag)
{
   // Physical constants
   const CONVERGE_precision_t R_univ = 8314.46;  // J/(kmol·K)
   const CONVERGE_precision_t PI_const = 3.14159265358979323846;
   
   // Calculate equilibrium mole fraction at droplet surface
   CONVERGE_precision_t chi_eq = P_sat / (P_gas + 1.0e-10);
   
   // Calculate Knudsen layer length scale
   // L_k = μ_g * sqrt(2π * R_univ * T_d / M_species) / (Sc * P_g)
   CONVERGE_precision_t sqrt_term = sqrt(2.0 * PI_const * R_univ * T_drop / M_species);
   CONVERGE_precision_t L_k = mu_gas * sqrt_term / (Sc_gas * P_gas + 1.0e-30);
   
   // Calculate non-dimensional correction parameter
   // ψ = 2 * L_k / d = L_k / R
   CONVERGE_precision_t psi = L_k / (radius + 1.0e-30);
   
   // Calculate evaporation rate parameter β
   // β = -(ρ_l * Pr * R / μ_g) * dR/dt
   CONVERGE_precision_t beta = -(rho_liquid * Pr_gas * radius / (mu_gas + 1.0e-30)) * drdt_prev;
   
   // Apply L-K correction to mole fraction
   // χ_neq = χ_eq / (1 + ψ·β)
   CONVERGE_precision_t chi_neq = chi_eq / (1.0 + psi * beta);
   
   // Clip to physical range [0, 0.9999]
   chi_neq = (chi_neq < 0.0) ? 0.0 : chi_neq;
   chi_neq = (chi_neq > 0.9999) ? 0.9999 : chi_neq;
   
   // Convert mole fraction to mass fraction
   // Y_v_s = (χ_neq * M_species) / (χ_neq * M_species + (1 - χ_neq) * M_gas)
   CONVERGE_precision_t numerator = chi_neq * M_species;
   CONVERGE_precision_t denominator = numerator + (1.0 - chi_neq) * M_gas_avg + 1.0e-30;
   CONVERGE_precision_t Y_v_s = numerator / denominator;
   
   // Clip to CONVERGE physical range
   Y_v_s = (Y_v_s < 1.0e-10) ? 1.0e-10 : Y_v_s;
   Y_v_s = (Y_v_s > 1.0) ? 1.0 : Y_v_s;
   
   // Diagnostic output to file (only on rank 0 to avoid conflicts)
   // Control with lk_diagnostic_flag in user_inputs.in (0=off, 1=on)
   // WARNING: Can create very large files (15+ GB) - use only when needed
   if(lk_diagnostic_flag == 1)
   {
      int rank;
      CONVERGE_mpi_comm_rank(&rank);
      
      if(rank == 0)  // Only rank 0 writes to file
      {
         // Write header on first call
         if(lk_diagnostic_header_written == 0)
         {
            FILE *fp_lk = fopen("lk_evap_diagnostics.txt", "w");
            if (fp_lk != NULL)
            {
               fprintf(fp_lk, "# Langmuir-Knudsen Evaporation Model Diagnostics\n");
               fprintf(fp_lk, "# ncyc = cycle number, time = simulation time [s]\n");
               fprintf(fp_lk, "# T_drop = droplet temperature [K], P_gas = gas pressure [Pa], P_sat = saturation pressure [Pa]\n");
               fprintf(fp_lk, "# R = radius [m], drdt = rate of radius change [m/s]\n");
               fprintf(fp_lk, "# chi_eq = equilibrium mole fraction, chi_neq = non-equilibrium mole fraction (LK corrected)\n");
               fprintf(fp_lk, "# L_k = Knudsen layer length [m], psi = L_k/R [-], beta = evaporation rate parameter [-]\n");
               fprintf(fp_lk, "# Y_v_s = surface vapor mass fraction (LK corrected), Y_inf = ambient vapor mass fraction\n");
               fprintf(fp_lk, "# Format: ncyc, time, T_drop, P_gas, P_sat, R, drdt, chi_eq, chi_neq, L_k, psi, beta, Y_v_s, Y_inf\n");
               fclose(fp_lk);
               lk_diagnostic_header_written = 1;
            }
         }
         
         // Append data
         FILE *fp_lk = fopen("lk_evap_diagnostics.txt", "a");
         if (fp_lk != NULL)
         {
            fprintf(fp_lk, "%ld, %.6e, %.3f, %.1f, %.1f, %.3e, %.3e, %.6f, %.6f, %.3e, %.6f, %.6f, %.6f, %.6f\n",
                    CONVERGE_ncyc(), CONVERGE_simulation_time(), T_drop, P_gas, P_sat,
                    radius, drdt_prev, chi_eq, chi_neq, L_k, psi, beta, Y_v_s, Y_inf);
            fclose(fp_lk);
         }
      }
   }
   
   return Y_v_s;
}

/**
 * @brief lk_residual_beta: computes the residual F(beta) = beta - (-drdt(beta) * rho_l * Pr_g * r / mu_g)
 */
static CONVERGE_precision_t lk_residual_beta(CONVERGE_precision_t beta_guess, void* ctx_ptr)
{
   LKCtx* ctx = (LKCtx*)ctx_ptr;
   
   // Deducible drdt from the trial beta
   // beta_guess = -(rho_l * Pr_g * R / mu_g) * drdt_guess
   CONVERGE_precision_t C_beta = (ctx->rho_liquid * ctx->pr_num * ctx->radius) / (ctx->mol_visc + 1.0e-30);
   CONVERGE_precision_t drdt_guess = -beta_guess / (C_beta + 1.0e-30);
   
   // 1. Calculate Y_1_star resulting from this trial drdt_guess
   CONVERGE_precision_t y1_star_trial = calculate_lk_y1_star(
      ctx->temp1,
      ctx->global_pressure,
      ctx->vapor_pres,
      ctx->radius,
      drdt_guess,
      ctx->mol_visc,
      ctx->rho_liquid,
      ctx->pr_num,
      ctx->sc_num,
      ctx->y1_inf,
      ctx->m_species,
      ctx->w_0,
      0 // No diagnostics inside residual
   );
   
   // 2. Re-evaluate Spalding mass transfer
   CONVERGE_precision_t bsub_d_trial = (y1_star_trial - ctx->y1_inf) / (1.0 - y1_star_trial + 1.0e-10);
   bsub_d_trial = (bsub_d_trial < 1.0e-15) ? 1.0e-15 : bsub_d_trial;
   
   // Evaluate drdt_eval
   // Based on standard formulation: drdt = -mass_trans_coeff * log(1 + B_M)
   // However mass_trans_coeff geometrically contains 1/(2 * r * rho), which cancels elegantly.
   // Utilizing mass_trans_coeff_geom = mass_trans_coeff * (2 * r * rho)
   CONVERGE_precision_t mass_trans_coeff = ctx->mass_trans_coeff_geom / (2.0 * ctx->radius * ctx->rho_liquid + 1.0e-30);
   CONVERGE_precision_t drdt_eval = -mass_trans_coeff * log(bsub_d_trial + 1.0);
   
   // Condensation safety
   if (drdt_eval >= 0.0) {
       drdt_eval = 0.0;
   }
   
   // Mass clipping considerations
   CONVERGE_precision_t evap_radius = ctx->radius + (drdt_eval * ctx->dt);
   if (evap_radius < 0.0) { evap_radius = 0.0; }
   
   CONVERGE_precision_t evap_mass_drop_1 = (4.0 / 3.0) * 3.14159265358979323846 * ctx->rho_liquid * (evap_radius*evap_radius*evap_radius);
   CONVERGE_precision_t dm_dt_eval = (evap_mass_drop_1 - ctx->mass_drop_tm1) / (ctx->dt + 1.0e-30);
   
   if ((-dm_dt_eval * ctx->dt) > ctx->evap_mass_drop_0) {
      dm_dt_eval = -ctx->evap_mass_drop_0 / (ctx->dt + 1.0e-30);
      // Back calculate drdt from dm_dt constraint
      CONVERGE_precision_t mass_new = ctx->mass_drop_tm1 + dm_dt_eval * ctx->dt;
      if (mass_new < 0.0) { mass_new = 0.0; }
      CONVERGE_precision_t r_new = pow((3.0 / 4.0) * mass_new / (3.14159265358979323846 * ctx->rho_liquid + 1.0e-30), 1.0/3.0);
      drdt_eval = (r_new - ctx->radius) / (ctx->dt + 1.0e-30);
   }
   
   // 3. Compare evaluated drdt against the one implicit to beta
   // F(beta) = beta - beta(drdt_eval)
   CONVERGE_precision_t beta_eval = -C_beta * drdt_eval;
   
   return beta_guess - beta_eval;
}
// /**
// * @brief fuller_diffusion_coef_nh3n2: Calculates the Fuller diffusion coefficient for NH3/N2
// * 
// * @param T: Gas Temperature in Kelvin 
// * @param P: Gas Pressure in Pa
// * 
// * @return Fuller diffusion coefficient in m^2/s
// */
// static CONVERGE_precision_t fuller_diffusion_coef_nh3n2(CONVERGE_precision_t T, CONVERGE_precision_t P)
// {
//    CONVERGE_precision_t fuller_coef = 0.00143;
//    CONVERGE_precision_t Vol_nh3 = 20.7;   //Properties of Gases and Liquids 5th ed. Table 11-1
//    CONVERGE_precision_t Vol_n2 = 18.5;   //Properties of Gases and Liquids 5th ed. Table 11-1
//    CONVERGE_precision_t M_NH3 = 17.031; //g/mol
//    CONVERGE_precision_t M_N2  = 28.014; //g/mol
//    CONVERGE_precision_t Vol_sum = CONVERGE_square(CONVERGE_cbrt(Vol_nh3) + CONVERGE_cbrt(Vol_n2));
//    CONVERGE_precision_t T_term = CONVERGE_pow(T, 1.75);
//    CONVERGE_precision_t P_term = P/1.0e5;
//    CONVERGE_precision_t M_ab = 1.0/(2 * ((1.0/M_NH3) + (1.0/M_N2)));
//    CONVERGE_precision_t M_ab_term = CONVERGE_sqrt(M_ab);


//    CONVERGE_precision_t D_fuller = fuller_coef * T_term / (P_term * Vol_sum * M_ab_term);

//    //Convert to m^2/s
//    D_fuller = D_fuller * 1.0e-4;
//    return D_fuller;
// }

/**
 * @brief solve_lk_bisection: Solves F(beta) = 0 via Bisection Method
 */
static CONVERGE_index_t solve_lk_bisection(CONVERGE_precision_t* beta_out,
                                           CONVERGE_precision_t beta_lo,
                                           CONVERGE_precision_t beta_hi,
                                           CONVERGE_precision_t tol_x,
                                           CONVERGE_precision_t tol_f,
                                           int max_iter,
                                           void* ctx)
{
   CONVERGE_precision_t f_lo = lk_residual_beta(beta_lo, ctx);
   CONVERGE_precision_t f_hi = lk_residual_beta(beta_hi, ctx);
   
   // Verify bracket
   int expansions = 0;
   while (f_lo * f_hi > 0.0 && expansions < 12) {
      beta_hi *= 10.0; // expand right
      f_hi = lk_residual_beta(beta_hi, ctx);
      expansions++;
   }
   
   if (f_lo * f_hi > 0.0) {
      // Bracket failed.
      return 1; 
   }
   
   CONVERGE_precision_t beta_mid, f_mid;
   for (int iter = 0; iter < max_iter; ++iter) {
      beta_mid = 0.5 * (beta_lo + beta_hi);
      f_mid = lk_residual_beta(beta_mid, ctx);
      
      if (fabs(beta_hi - beta_lo) < tol_x || fabs(f_mid) < tol_f) {
         *beta_out = beta_mid;
         return 0; // Success
      }
      
      if (f_lo * f_mid < 0.0) {
         beta_hi = beta_mid;
         f_hi = f_mid;
      } else {
         beta_lo = beta_mid;
         f_lo = f_mid;
      }
   }
   
   // Failed to converge within max_iter
   *beta_out = beta_mid;
   return 2; 
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
