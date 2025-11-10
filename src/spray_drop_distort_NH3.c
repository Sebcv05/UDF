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
#include <spray_break.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<unistd.h>
#include <assert.h>
#include <complex.h>
#include <time.h>
#include <CubicSolver.h>
#include <BreakupCriterion.h>
#include <TABDistort.h>
#include <Breakup.h>
#include <BubbleDensityNH3.h>
#include <PsatNH3.h>
#include <TsatNH3.h>
#include <DGRE_NH3.h>
#include <Geometry.h>
//#include <mpi.h>
//#include <DestroyTables.h>
#include <Vb.h>
#include <RPE_euler.h>
#include <globals.h>


/// @brief UDF Thermal Breakup Model with all thermal properties for ammonia
/// @param mesh 
/// @param cloud 
/// @param spray_cloud_list 
/// @param i_pc 
/// @param node_index 
/// @param global_density 
/// @param global_viscosity 
/// @param parcel_counter 
/// @param sp 
static void spray_distort_cell_NH3(CONVERGE_mesh_t mesh, CONVERGE_cloud_t cloud, CONVERGE_cloud_list_t spray_cloud_list, CONVERGE_index_t i_pc, CONVERGE_index_t node_index, const CONVERGE_precision_t *global_density, const CONVERGE_precision_t *global_viscosity, CONVERGE_index_t parcel_counter,CONVERGE_species_t sp);
static void sync_child_velocity(CONVERGE_cloud_t cloud, CONVERGE_cloud_list_t spray_cloud_list, CONVERGE_index_t i_pc, CONVERGE_index_t node_index);
static void init_tables(CONVERGE_species_t species);
static void destroy_tables(CONVERGE_species_t species);

// static void init_tables(CONVERGE_species_t species);
// static void destroy_tables(CONVERGE_species_t species);
static CONVERGE_table_t *pvap_table = NULL;
static CONVERGE_table_t *visc_table = NULL;
static CONVERGE_table_t *cond_table = NULL;
static CONVERGE_table_t *cp_table   = NULL;
static CONVERGE_table_t *rho_table  = NULL;
static CONVERGE_table_t *h_table    = NULL;
static CONVERGE_table_t *evap_species_h_table;
static CONVERGE_table_t *evap_species_sensible_h_table;

// Profiling accumulators
static double prof_geom=0.0, prof_dgre=0.0, prof_bc=0.0, prof_break=0.0, prof_bubble=0.0;
static int last_cycle = -1;
static int spray_params_logged = 0;
/**********************************************************************/
/*                                                                    */
/* Name: user_drop_distort, user_drop_distort_cell                    */
/*                                                                    */
/* Description: user_drop_distort includes a drop distortion          */
/* calculation based on the TAB model from: The TAB Method for        */
/* Numerical Calculation of Spray Droplet Breakup, O'Rourke and       */
/* Amsden SAE 872089.                                                 */
/*                                                                    */
/* Called when: user_distort_flag=1 in udf.in                         */
/*                                                                    */
/* Variables passed in: none                                          */
/*                                                                    */
/* Variables passed back: none                                        */
/*                                                                    */
/* Subroutines called: BubbleDensity.c, Psat.c, Tsat.c  */
/*     Geometry.c, DGRE.c, BreakupCriterion.c, Breakaup.c             */
/**********************************************************************/
static CONVERGE_table_t *hvap_table = NULL;



CONVERGE_UDF(drop_distort, IN(FIELD(CONVERGE_precision_t *, density), VALUE(CONVERGE_mesh_t, mesh), FIELD(CONVERGE_precision_t *, pressure), FIELD(CONVERGE_precision_t *, temperature), FIELD(CONVERGE_precision_t *, gas_mol_viscosity)), OUT(CONVERGE_VOID))
{

   CONVERGE_precision_t distort_start = CONVERGE_mpi_wtime();

   CONVERGE_index_t parcel_counter = 0;
   const CONVERGE_cloud_list_t spray_cloud_list = CONVERGE_mesh_get_spray_cloud_list(mesh);
   global_pressure = pressure;
   global_temperature = temperature;
   global_density = density;
   global_mol_viscosity = gas_mol_viscosity;
   CONVERGE_cloud_t cloud;
   CONVERGE_iterator_t cl_it;
   CONVERGE_species_t sp = CONVERGE_mesh_species(mesh);
   CONVERGE_cloud_list_iterator_create(spray_cloud_list, &cl_it);
   init_tables(sp);
   for (CONVERGE_index_t i_pc = CONVERGE_iterator_first(cl_it); i_pc != -1; i_pc = CONVERGE_iterator_next(cl_it))
   {
      cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, i_pc);
      const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);
      // printf("\ninitializing tables");
     
      // printf("\n spray_distort_cell...");
    
   
      //Timing 
      
      spray_distort_cell_NH3(mesh, cloud, spray_cloud_list, i_pc, node_index, density, gas_mol_viscosity, parcel_counter,sp);
     
     // printf("\ndiff = %es",sdc_diff);
    //  printf("\n destroying tables...");
      

   }
   destroy_tables(sp);
  
   //Get rank 
   // CONVERGE_int_t rank;
   //  CONVERGE_mpi_comm_rank(&rank);
   // if(rank==0)
   // {
   // printf("cloud loop time = %f\n",cl_diff);
   // }
   CONVERGE_iterator_destroy(&cl_it);

   CONVERGE_precision_t end_time = CONVERGE_mpi_wtime();
   // CONVERGE_precision_t total_time = end_time - start_time;

     // printf("\nTotal time: %f\n",total_time);
   
     CONVERGE_precision_t distort_end = CONVERGE_mpi_wtime();
     CONVERGE_int_t rank;
     CONVERGE_mpi_comm_rank(&rank);
     CONVERGE_precision_t distort_diff = distort_end - distort_start;
   //   printf("\n spray distort rank %d runtime = %.2e (s)\n",rank,distort_diff);
}



static void spray_distort_cell_NH3(CONVERGE_mesh_t mesh, CONVERGE_cloud_t cloud, CONVERGE_cloud_list_t spray_cloud_list, CONVERGE_index_t i_pc, CONVERGE_index_t node_index, const CONVERGE_precision_t *global_density, const CONVERGE_precision_t *global_viscosity, CONVERGE_index_t parcel_counter,CONVERGE_species_t sp)
{
   /*  fprintf(stderr,"mesh size = %i\n",sizeof(mesh));
            fprintf(stderr,"cloud size = %i\n",sizeof(cloud));
                  fprintf(stderr,"mesh size = %i\n",sizeof(spray_cloud_list));
*/



   //MB chck
   CONVERGE_precision_t mass_before, mass_after;
   
   // Breakup event counter
   static int total_breakups = 0;
   static int last_reported_cycle = -1;
   int breakups_this_call = 0;
   
   // Single parcel tracking for detailed diagnostics
   static FILE* parcel_track_file = NULL;
   static int tracking_initialized = 0;
   static int tracked_parcel_id = -1;  // Will be set to first parent parcel we encounter



   // printf("\n0");
   CONVERGE_size_t num_parcels_in_cloud = CONVERGE_cloud_size(cloud);
   struct ParcelCloud old_parcel_cloud, new_parcel_cloud;

   load_user_cloud(&old_parcel_cloud, cloud);

   CONVERGE_precision_t pre_pl = CONVERGE_mpi_wtime();
   CONVERGE_precision_t pre_TAB,post_TAB,pre_DGRE,post_DGRE,pre_Geom,post_Geom,pre_break,post_break,pre_bc,post_bc,pre_pbr,post_pbr,sopl,eopl;
   // printf("\n 0.1");
   // printf("starting loop over parcels in cloud\n");
   mass_before = 0;
   mass_after = 0;
   for (int p_idx = 0; p_idx < num_parcels_in_cloud; p_idx++)
   {
      // if(CONVERGE_simulation_time_sec() > 1.0e-4)
      // {
      //Timing 
         pre_TAB = 0.0; post_TAB = 0.0; pre_DGRE=0.0;post_DGRE=0.0;pre_Geom=0.0;post_Geom=0.0;pre_break=0.0;post_break=0.0;pre_bc=0.0;pre_bc=0.0;pre_pbr=0.0;post_bc=0.0;
         sopl = CONVERGE_mpi_wtime();

         mass_before= mass_before + (1.33333 * PI * old_parcel_cloud.num_drop[p_idx]*CONVERGE_cube(old_parcel_cloud.radius[p_idx]));
         old_parcel_cloud.m0[p_idx] = (1.33333 * PI * CONVERGE_cube(old_parcel_cloud.radius[p_idx])*old_parcel_cloud.num_drop[p_idx]);
        // printf("\n before num_drop = %e rad = %e",old_parcel_cloud.num_drop[p_idx],old_parcel_cloud.radius[p_idx]);

      CONVERGE_precision_t time_start = CONVERGE_simulation_time_sec();
      // printf("\n time_start = %e",time_start);
// printf("\n0.2");
   CONVERGE_size_t num_parcels_before = CONVERGE_cloud_list_num_parcels(spray_cloud_list);
   if (num_parcels_before <= 0)
   {
      // printf("OUT: \n");
      return;
   }

   // printf("\n0.3");
   CONVERGE_size_t num_parcel_species = CONVERGE_species_num_parcel(sp);

//   printf("\n1"); 

   CONVERGE_index_t new_pc_idx, new_p_idx;
   // printf("user cloud loaded\n");
   
   // printf("loaded global variables\n");
   // printf("rho_v = %e\n",rho_v);
   // printf("mu_v = %e \n",mu_v);
   //  Simulation Metadata
   // printf("mesh variables loaded\n");
   CONVERGE_precision_t dt = CONVERGE_simulation_dt();
   // Get high-level API containers
   num_gas_species = CONVERGE_species_num_gas(sp);
   num_parcel_species = CONVERGE_species_num_parcel(sp);

   // Old table lookup vars
   CONVERGE_precision_t average_hvap = 0.0;
   CONVERGE_precision_t hvap;

   // Initialize variables and tables.
   


   CONVERGE_int_t theskyisblue = 1;          // it is 
   CONVERGE_int_t theskyisgreen = 0;         //it isn't, other than in Skinner's Kitchen
   /*--------------------------------------------------------------*/


   // printf("iterator created\n");
   //  Local variables
   CONVERGE_precision_t mu;       // Liquid Viscosity
   CONVERGE_precision_t sigma;    // Surface Tension
   CONVERGE_precision_t Td;       // Droplet Temperature
   CONVERGE_precision_t Rb;       // Bubble Radius
   CONVERGE_precision_t Rb_0;     // Previous Bubble Radius
   CONVERGE_precision_t Rb_temp;  // Temporary holder for Rb
   CONVERGE_precision_t Vb_tm1;   // Previous Bubble Velocity
   CONVERGE_precision_t rho_l;    // Droplet Density
   CONVERGE_precision_t k;        // Suface Viscosity
   CONVERGE_precision_t H;        // Enthalpy of Vaporisation
   CONVERGE_precision_t A;        // Constant for growth = sqrt(2dP/3)i
   CONVERGE_precision_t t_parcel; // Parcel Lifetime
   CONVERGE_precision_t csubp_l;  // Liquid specific heat
   CONVERGE_precision_t alpha_l;  // Constant from Plesset & Zwick
   CONVERGE_precision_t tau;      // Non-dimensional Time
   CONVERGE_precision_t NDR;      // Non dimensional Radius
   CONVERGE_precision_t rho_b;    // Bubble density
   CONVERGE_precision_t tab_om;
   // printf("initiated local variables\n");
   //  Mesh Vars

   CONVERGE_precision_t P_amb = global_pressure[node_index];
   CONVERGE_precision_t T_amb = global_temperature[node_index];
   CONVERGE_precision_t csubp = global_csubp[node_index];
   CONVERGE_precision_t rho_v = global_density[node_index];
   CONVERGE_precision_t mu_v = global_mol_viscosity[node_index] * rho_v;

      // old_parcel_cloud.tbreak_rt[p_idx] = old_parcel_cloud.temp[p_idx];
      parcel_counter++;
      // Store thermal breakup flag in from_nozzle so it can be exported in .h5s
      //
      // old_parcel_cloud.from_nozzle[p_idx] = old_parcel_cloud.thermal_breakup_flag[p_idx];
      // printf("p_idx = %i\n",p_idx);
      //   Populate local variables

      sigma = old_parcel_cloud.surf_ten[p_idx];
      Td = old_parcel_cloud.temp[p_idx];
      Rb = old_parcel_cloud.r_bubble[p_idx];
      Rb_0 = old_parcel_cloud.r_bubble_0[p_idx];
      Vb_tm1 = old_parcel_cloud.v_bubble_tm1[p_idx];
      k = 0; // Initial Assumption
      t_parcel = old_parcel_cloud.lifetime[p_idx];
      CONVERGE_precision_t g_den = rho_v;
      CONVERGE_precision_t g_pressure = P_amb;
      
      // DIAGNOSTIC: Check radius at very start of loop
      static int radius_diag_count = 0;
      if (radius_diag_count < 5 && old_parcel_cloud.is_child[p_idx] == 0) {
         printf("[DISTORT_START] p_idx=%d, radius=%.6e m, lifetime=%.6e s, Td=%.2f K, is_child=%d\n",
                p_idx, old_parcel_cloud.radius[p_idx], old_parcel_cloud.lifetime[p_idx], 
                Td, old_parcel_cloud.is_child[p_idx]);
         radius_diag_count++;
      }
      
      // Initialize single parcel tracking
      if (!tracking_initialized && old_parcel_cloud.is_child[p_idx] == 0 && old_parcel_cloud.pbt[p_idx] == 1) {
         tracked_parcel_id = p_idx;
         tracking_initialized = 1;
         parcel_track_file = fopen("tracked_parcel.csv", "w");
         if (parcel_track_file) {
            fprintf(parcel_track_file, "time,lifetime,p_idx,R_drop,R_bubble,Rdot,T_drop,Pb,P_amb,thermal_breakup_flag,kb\n");
            printf("[TRACKING] Started tracking parcel %d\n", tracked_parcel_id);
         }
      }
      
      // Note: Detailed logging happens later after kb is calculated
      
      // Calculate Saturation Pressure from Antoine's Equation
      CONVERGE_precision_t P_sat;
      if(Td > 300.0){
         printf("\n\n before P_sat, Td = %f, radius = %e, p_idx = %li, tbt = %i, pbt = %i, thermal_breakup_flag = %d",Td,old_parcel_cloud.radius[p_idx],old_parcel_cloud.tbt[p_idx],old_parcel_cloud.pbt[p_idx],old_parcel_cloud.thermal_breakup_flag[p_idx]);
         printf("\nRemoving Parcel\n");
         old_parcel_cloud.num_drop[p_idx] = 0.0;
         old_parcel_cloud.radius[p_idx] = 0.0;
         old_parcel_cloud.r_bubble[p_idx] = 0;
         old_parcel_cloud.r_bubble_0[p_idx] = 0;
         old_parcel_cloud.v_bubble_tm1[p_idx] = 0;
         old_parcel_cloud.lifetime[p_idx] = 0;
         old_parcel_cloud.tbt[p_idx] = 0;
         old_parcel_cloud.pbt[p_idx] = 0;
         old_parcel_cloud.thermal_breakup_flag[p_idx] = 0;
         old_parcel_cloud.int_omega[p_idx] = 0;
         old_parcel_cloud.int_omega[p_idx] = 0;
         old_parcel_cloud.m0[p_idx] = 0;
         old_parcel_cloud.dgre_cycle_count[p_idx] = 0;
         old_parcel_cloud.r_bubble_0[p_idx] = old_parcel_cloud.r_bubble[p_idx];
         old_parcel_cloud.r_therm[p_idx] = old_parcel_cloud.radius[p_idx];
         continue;
      }
      Saturation_PressureNH3(Td, &P_sat);
      // Unit tests for saturation pressure function
      CONVERGE_precision_t Psattest1, Psattest2, Psattest3, Ttest1, Ttest2, Ttest3;
      Ttest1 = 300.274;

    Saturation_PressureNH3(Ttest1, &Psattest1);
      if (round(Psattest1) != 1061549)
      {
         printf("\nP_sat function failed unit test, aborting...");
         printf("\n %f != 77188", floor(Psattest1 * 1e5) / 1e5);
         CONVERGE_mpi_abort();
      }
      // Unit tests 2 & 3 -- only uncomment to check program aborts when T outside of range
      //    Ttest2=280;
      //    Ttest3=400;
      //   // Saturation_PressureNH3(Ttest2,&Psattest2);
      //    Saturation_PressureNH3(Ttest3,&Psattest3);

      if (P_sat < 1)
      {
         printf("\nParcel Temperature outside of range for P_sat Correlation, continuing...");
         old_parcel_cloud.v_bubble[p_idx] =0.0;
         old_parcel_cloud.omega[p_idx] = 0.0;
         old_parcel_cloud.eta_drop[p_idx] = 0.0;
         old_parcel_cloud.int_omega[p_idx] =0.0;
         old_parcel_cloud.int_omega[p_idx] = 0.0;
         old_parcel_cloud.m0[p_idx] = (1.33333 * PI * CONVERGE_cube(old_parcel_cloud.radius[p_idx])*old_parcel_cloud.num_drop[p_idx]);
         old_parcel_cloud.dgre_cycle_count[p_idx] = 0;
         old_parcel_cloud.r_bubble_0[p_idx] = old_parcel_cloud.r_bubble[p_idx];
         old_parcel_cloud.r_therm[p_idx] = old_parcel_cloud.radius[p_idx];
         continue;

      }
      rho_b = bubble_densityNH3(P_sat, Td); // Not sure about using this to estimate rho_b
    

      
         CONVERGE_precision_t P_sat_new;
         Saturation_PressureNH3(old_parcel_cloud.temp[p_idx],&P_sat_new);
         // old_parcel_cloud.r_bubble[p_idx]= 2.0 * old_parcel_cloud.surf_ten[p_idx] / (P_sat_new - global_pressure[node_index]);
         if(old_parcel_cloud.r_bubble[p_idx]<0.0)     //not superhearted anymore, deactivate thermal model for this parcel
         {
            old_parcel_cloud.r_bubble[p_idx]=0.0;
            old_parcel_cloud.thermal_breakup_flag[p_idx]=999;
            old_parcel_cloud.pbt[p_idx]=0;
         
         }


   


      theskyisblue = 1;    //it is 
      theskyisgreen = 0;   // it is not
      if (theskyisblue)
      {

         // Pre-breakup routine
         pre_pbr = CONVERGE_mpi_wtime();

    
      

         if (old_parcel_cloud.thermal_breakup_flag[p_idx] <0 && old_parcel_cloud.pbt[p_idx]==1)
         {
            if(old_parcel_cloud.thermal_breakup_flag[p_idx]>0)
            {
                           printf("\n tbf at start of loop is %i",old_parcel_cloud.thermal_breakup_flag[p_idx]);

            }
            if(old_parcel_cloud.radius[p_idx]> old_parcel_cloud.r_drop_0[p_idx]*1.5)
            {
               old_parcel_cloud.thermal_breakup_flag[p_idx] = 1;
               old_parcel_cloud.tbt[p_idx] = 1;
               old_parcel_cloud.pbt[p_idx] = 0;
               old_parcel_cloud.thermal_breakup_flag[p_idx] = 6;
               continue;
               // printf("\n tbf at start of loop is %i",old_parcel_cloud.thermal_breakup_flag[p_idx]);
            
               // printf("\n tbf at start of loop is %i",old_parcel_cloud.thermal_breakup_flag[p_idx]);
            }

            //Calculate Species dependent properites
            CONVERGE_precision_t average_hvap = 0.0;
            csubp_l = 0.0;
            for (CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
            {
               average_hvap += old_parcel_cloud.mfrac[p_idx * num_parcel_species + isp] *  (CONVERGE_table_lookup(hvap_table[isp], Td) + 1.0e-3);
              csubp_l +=old_parcel_cloud.mfrac[p_idx * num_parcel_species + isp] * CONVERGE_table_lookup(cp_table[isp], Td);
            }
            // NaN avoider
            if ((1 + csubp_l * (T_amb - Td) / average_hvap) <= 0.0)
            {
            }
            H = average_hvap;

            /*    Use to estimate Ja number
            CONVERGE_precision_t T_saturation = T_sat(global_pressure[node_index]);
            CONVERGE_precision_t dT = old_parcel_cloud.temp[p_idx] - T_saturation;
            CONVERGE_precision_t Ja = dT * csubp_l/H; 
            //printf("\nH = %e cp = %e cp/H = %e dT = %e Ja = %f",H,csubp_l,csubp_l/H,dT,Ja);
            // printf("\n H = %e	T = %e	csubp_l = %e",average_hvap,Td,csubp_l);
               */

            // --- Synchronize bubble and droplet size in case of KH-RT or other breakup shrinkage ---
            if (old_parcel_cloud.r_bubble[p_idx] > old_parcel_cloud.radius[p_idx]) {
                CONVERGE_precision_t ratio = old_parcel_cloud.r_bubble_0[p_idx] / old_parcel_cloud.r_drop_0[p_idx];
                old_parcel_cloud.r_bubble[p_idx] = ratio * old_parcel_cloud.radius[p_idx];
          
                old_parcel_cloud.r_bubble_0[p_idx] = old_parcel_cloud.r_bubble[p_idx];
               //  printf("\n[Bubble sync] KH-RT or secondary breakup reduced r_drop; rescaled r_bubble to maintain ratio.\n");
            }
            //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // Start of sub cycle loop
            //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
              // target physical substep size (s)
              // Reduced from 1e-10 to 1e-12 during debugging - now safe_divide bug is fixed, can restore
               const CONVERGE_precision_t dt_sub_target = 1.0e-9;

               // current CFD timestep
               CONVERGE_precision_t dt_global = CONVERGE_simulation_dt();

               // compute how many substeps are needed
               int n_sub = (int)ceil(dt_global / dt_sub_target);
               if (n_sub < 1) n_sub = 1;

               // recompute actual substep (in case dt_global not exact multiple)
               CONVERGE_precision_t dt_sub = dt_global / (CONVERGE_precision_t)n_sub;
               for (int sub = 0; sub < n_sub; ++sub) {

               //RPE Euler solver - updates v_bubble, r_bubble, temp
               CONVERGE_precision_t t0 = CONVERGE_mpi_wtime();
               RPE_euler_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub, 
                                hvap_table, cp_table, num_parcel_species);
               prof_bubble += CONVERGE_mpi_wtime() - t0;
               
               // Check if bubble growth stopped
               if(old_parcel_cloud.v_bubble[p_idx]<1.0e-10)
               {
                  old_parcel_cloud.pbt[p_idx] = 0;
                  old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;
                  break;
               }

               //Load Rb
               Rb = old_parcel_cloud.r_bubble[p_idx];
               Rb_temp = Rb;

            // RPE_euler already updated r_bubble, just do safety checks
            Rb = old_parcel_cloud.r_bubble[p_idx];
            
            if (isnan(Rb))
            {
               printf("\nRb is NaN after RPE solver\n");
               CONVERGE_mpi_abort();
            }
            
            if (Rb < 0.0)
            {
               Rb = 0.0;
               old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;
               old_parcel_cloud.r_bubble[p_idx] = Rb;
               printf("Rb negative after RPE solver\n");
               CONVERGE_mpi_abort();
            }
            
            if (Rb > old_parcel_cloud.radius[p_idx])
               {
                  old_parcel_cloud.thermal_breakup_flag[p_idx] = 5;
                  old_parcel_cloud.tbt[p_idx]=1;
                  old_parcel_cloud.r_bubble[p_idx] = 0.8 * old_parcel_cloud.radius[p_idx];
                  break;
               }
               else
               {
                  //Save Rb 
                  old_parcel_cloud.r_bubble[p_idx] = Rb;


            // Update droplet radius
            t0 = CONVERGE_mpi_wtime();
            CONVERGE_precision_t rdrop_before_geometry = old_parcel_cloud.radius[p_idx];

            Geometry(&old_parcel_cloud, p_idx, dt_sub);
            prof_geom += CONVERGE_mpi_wtime() - t0;
            CONVERGE_precision_t rdrop_after_geometry = old_parcel_cloud.radius[p_idx];
                  if(rdrop_after_geometry< rdrop_before_geometry){

                     printf("\nGeometry has shrunk parcel, rdrop_before_geometry = %e, rdrop_after_geometry = %e\n", rdrop_before_geometry, rdrop_after_geometry);
                     printf("\nAborting...");
                     CONVERGE_mpi_abort();
                  }


            CONVERGE_precision_t g_den = global_density[node_index];
      
   
            t0 = CONVERGE_mpi_wtime();
            if(old_parcel_cloud.r_bubble[p_idx] > 0.02 * old_parcel_cloud.radius[p_idx]){
            DGRE_NH3(&old_parcel_cloud, p_idx, g_den);
            }
            prof_dgre += CONVERGE_mpi_wtime() - t0;
            // Update Breakup Criteria
            //  printf("calling BreakupCriterion\n");
            t0 = CONVERGE_mpi_wtime();
            CONVERGE_precision_t kb = BreakupCriterion(&old_parcel_cloud, p_idx, dt_sub);
            if(!spray_params_logged)
            {
               spray_params_logged = 1;
            }
            prof_bc += CONVERGE_mpi_wtime() - t0;
            
            // Log kb for tracked parcel
            if (tracking_initialized && p_idx == tracked_parcel_id && parcel_track_file) {
               CONVERGE_precision_t sim_time = CONVERGE_simulation_time_sec();
               CONVERGE_precision_t Pb_track;
               Saturation_PressureNH3(Td, &Pb_track);
               fprintf(parcel_track_file, "%.12e,%.12e,%d,%.12e,%.12e,%.12e,%.6f,%.6e,%.6e,%d,%.6e\n",
                       sim_time,
                       t_parcel,
                       p_idx,
                       old_parcel_cloud.radius[p_idx],
                       Rb,
                       old_parcel_cloud.v_bubble[p_idx],
                       Td,
                       Pb_track,
                       P_amb,
                       old_parcel_cloud.thermal_breakup_flag[p_idx],
                       kb);
               fflush(parcel_track_file);
            }
         

   
   
            if (kb > kb_threshold)
            {
               // printf("\n Breakup happening due to kb > 1.0, kb = %e, rb = %e, r_drop = %e, vb = %e\n",kb,Rb,old_parcel_cloud.radius[p_idx],old_parcel_cloud.v_bubble[p_idx]);
               old_parcel_cloud.thermal_breakup_flag[p_idx] = 3;
               old_parcel_cloud.tbt[p_idx] = 1;
               old_parcel_cloud.r_bubble[p_idx] = Rb;
               old_parcel_cloud.eta_drop[p_idx] = 0;
               break;
            }
   
            // Store eta
            // old_parcel_cloud.eta_drop[p_idx] += eta; // Update eta
   
               //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
               // End of sub cycle loop
               //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
   
   
         
         

            } // Delimiter for if(r_bubble>r_drop)
         } // Delimiter for sub cycle loop


         } // Delimiter for if thermal_breakup_flag=0
         //********************** BREAKUP ROUTINE ***************************************//
         // old_parcel_cloud.from_nozzle[p_idx] = old_parcel_cloud.thermal_breakup_flag[p_idx];
         post_pbr = CONVERGE_mpi_wtime();
         if (old_parcel_cloud.tbt[p_idx] && old_parcel_cloud.thermal_breakup_flag[p_idx]!=4)
         {
              breakups_this_call++;  // Count breakup event
              
              // Calculate bubble pressure at breakup
              CONVERGE_precision_t R_final = old_parcel_cloud.r_bubble[p_idx];
              CONVERGE_precision_t T_final = old_parcel_cloud.temp[p_idx];
              CONVERGE_precision_t Vb_final = (4.0/3.0) * PI * R_final * R_final * R_final;
              
              // Get bubble pressure from saturation at droplet temperature
              CONVERGE_precision_t Pb_final;
              Saturation_PressureNH3(T_final, &Pb_final);
              
              // Detailed breakup diagnostic
              CONVERGE_precision_t Vb_final = (4.0/3.0) * PI * R_final * R_final * R_final;
              CONVERGE_precision_t m_b_final = old_parcel_cloud.m_bubble[p_idx];
              CONVERGE_precision_t rho_v_final = (Vb_final > 1e-30) ? (m_b_final / Vb_final) : 0.0;
              CONVERGE_precision_t Pb_actual = (rho_v_final > 1e-6) ? (rho_v_final * 488.2 * T_final) : 0.0;
              
              printf("[BREAKUP] Parcel %d: lifetime=%.6e s, Rdot=%.6e m/s, T=%.2f K, Pb_eq=%.3e Pa, Pb_actual=%.3e Pa, R_bubble=%.6e m, m_b=%.3e kg, rho_v=%.3e kg/m3\n",
                     p_idx,
                     old_parcel_cloud.lifetime[p_idx],
                     old_parcel_cloud.v_bubble[p_idx],
                     T_final,
                     Pb_final,
                     Pb_actual,
                     R_final,
                     m_b_final,
                     rho_v_final);
              
              // If this is the tracked parcel, close the tracking file
              if (tracking_initialized && p_idx == tracked_parcel_id && parcel_track_file) {
                  fprintf(parcel_track_file, "# BREAKUP OCCURRED AT THIS TIMESTEP\n");
                  fclose(parcel_track_file);
                  parcel_track_file = NULL;
                  printf("[TRACKING] Breakup detected for tracked parcel %d - tracking file closed\n", tracked_parcel_id);
              }
              
              // printf("\n breakup tbf = %i",old_parcel_cloud.thermal_breakup_flag[p_idx]);
               
               pre_break = CONVERGE_mpi_wtime();
               // Check cloud size before breakup
               CONVERGE_index_t cloud_size_before_break = CONVERGE_cloud_size(cloud);
               // printf("\nBreakup.cloud_size_before_break = %i", cloud_size_before_break);
              CONVERGE_precision_t t0 = CONVERGE_mpi_wtime();
               Breakup(&old_parcel_cloud, p_idx, cloud);
               prof_break += CONVERGE_mpi_wtime() - t0;
               // printf("\nBreakup.cloud_size_after_break = %i\n\n\n\n", CONVERGE_cloud_size(cloud));
               post_break = CONVERGE_mpi_wtime();
         }
  
        
  
      } // End of if(theskyisblue)
             
            eopl = CONVERGE_mpi_wtime();

         CONVERGE_precision_t post_pl = CONVERGE_mpi_wtime();
         CONVERGE_precision_t pl_diff = eopl - sopl;
         CONVERGE_precision_t tab_diff = post_TAB - pre_TAB;
         CONVERGE_precision_t tab_frac = 100.00* tab_diff / pl_diff;
         CONVERGE_precision_t geom_diff = post_Geom - pre_Geom;
         CONVERGE_precision_t dgre_diff = post_DGRE - pre_DGRE;
         CONVERGE_precision_t geom_frac, dgre_frac,break_diff,break_frac,bc_diff,bc_frac,pbr_diff,pbr_frac,init_diff,init_frac,bub_diff,bub_frac;
         break_diff = post_break - pre_break;
         if(break_diff<0.0)
         {
            break_diff = 0.0;
         }
         break_frac = 100.00 * break_diff/ pl_diff;

         bc_diff = post_bc - pre_bc;
         
         pbr_diff = post_pbr - pre_pbr;
         pbr_frac = 100.00 * pbr_diff / pl_diff; 
         init_diff = pre_TAB - sopl;
         init_frac = 100.00 * init_diff /pl_diff;
         bub_diff = pre_Geom - pre_pbr;
         bub_frac = 100.00 * bub_diff / pbr_diff;
         bc_frac = 100.00 * bc_diff / pbr_diff; 
         geom_frac = 100.00 * geom_diff / pbr_diff;
         dgre_frac = 100.00 * dgre_diff / pbr_diff;
      //  if(old_parcel_cloud.thermal_breakup_flag[p_idx] ==4 && old_parcel_cloud.tbt[p_idx]==0)
      //  {
      //    printf("\npl_time = %f",pl_diff);
      //    printf("\n     dt_init    = %e  frac = %f\%",init_diff,init_frac);
      //    printf("\n     dt_TAB     = %e  frac= %f\%",tab_diff,tab_frac);
      //   // if(old_parcel_cloud.thermal_breakup_flag[p_idx]<0)
      //    {
      //    printf("\n     dt_pbr     = %e frac= %f\%",pbr_diff,pbr_frac);
      //    printf("\n           dt_bubble  = %e frac= %f\%",bub_diff,bub_frac);
      //    printf("\n           dt_geom    = %e frac= %f\%",geom_diff,geom_frac);
      //    printf("\n           dt_dgre    = %e frac= %f\%",dgre_diff,dgre_frac);
      //    printf("\n           dt_bc      = %e frac= %f\%",bc_diff,bc_frac);
      //    }
      //    printf("\n     dt_breakup = %e frac= %f\%",break_diff,break_frac);
      //  }
    mass_after= mass_after + (1.33333 * PI * old_parcel_cloud.num_drop[p_idx]*CONVERGE_cube(old_parcel_cloud.radius[p_idx]));       
     // printf("\n after num_drop = %e rad = %e",old_parcel_cloud.num_drop[p_idx],old_parcel_cloud.radius[p_idx]);
      // } //time limieter >0.1ms 
   }    // End of parcel loop
   
   // Report breakup events
   int ncyc = CONVERGE_ncyc();
   total_breakups += breakups_this_call;
   
   if (ncyc != last_reported_cycle && breakups_this_call > 0) {
      int rank;
      CONVERGE_mpi_comm_rank(&rank);
      CONVERGE_precision_t sim_time = CONVERGE_simulation_time();
      printf("[BREAKUP] Cycle %d, Rank %d, Time %.6e s: %d breakups this cloud (Total: %d)\n",
             ncyc, rank, sim_time, breakups_this_call, total_breakups);
      last_reported_cycle = ncyc;
   }
   
   if (ncyc != last_cycle) {
    int rank;
    CONVERGE_mpi_comm_rank(&rank);
  
   //     // Print every rank’s data
   //  printf("Rank %d, Cycle %d profiling (s): Bubble=%f, Geometry=%f, DGRE=%f, BreakupCriterion=%f, Breakup=%f\n",
   //    rank, ncyc, prof_bubble, prof_geom, prof_dgre, prof_bc, prof_break);
   //    CONVERGE_precision_t prof_total = prof_bubble+prof_geom+prof_dgre+prof_bc+prof_break;
   //  // Print percentages
   //  printf("Rank %d, Cycle %d profiling (%%): Bubble=%f, Geometry=%f, DGRE=%f, BreakupCriterion=%f, Breakup=%f\n",
   //    rank, ncyc, 100.0f * prof_bubble / prof_total, 100.0f * prof_geom / prof_total, 100.0f * prof_dgre / prof_total, 100.0f * prof_bc / prof_total, 100.0f * prof_break / prof_total);



    // reset accumulators
    prof_geom = prof_dgre = prof_bc = prof_break = prof_bubble = 0.0;
    last_cycle = ncyc;
}
    

    CONVERGE_precision_t end_time = CONVERGE_mpi_wtime();
   //  CONVERGE_precision_t total_time = end_time - start_time;
   //  printf("\nTotal time: %f\n",total_time);
   
   // int rank;
   // CONVERGE_mpi_comm_rank(&rank);

   // printf("end of spray break cell\n");




   if(mass_after>mass_before)
   {
      CONVERGE_precision_t increase = mass_after-mass_before;
      if(mass_after/mass_before> 1.01)
      {
         printf("\n mass_before = %e mass_after = %e, increase = %e",mass_before,mass_after,mass_after-mass_before);

      }

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

   // evap_species_h_table =
   //    CONVERGE_get_species_prop_table(TEMPERATURE_TABLE_ENTHALPY_ID, CONVERGE_spray_evap_species_index());
   // evap_species_sensible_h_table =
   //    CONVERGE_get_species_prop_table(TEMPERATURE_TABLE_SENSIBLE_ENTHALPY_ID, CONVERGE_spray_evap_species_index());

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
   //CONVERGE_table_destroy(&evap_species_h_table);
   //CONVERGE_table_destroy(&evap_species_sensible_h_table);

   // Destroy local gas tables
   unload_species_tables(gas_species_it, &evap_species_h_table);
   unload_species_tables(gas_species_it, &evap_species_sensible_h_table);
   unload_species_tables(gas_species_it, &visc_table);
   unload_species_tables(gas_species_it, &cond_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
   CONVERGE_iterator_destroy(&gas_species_it);
}



CONVERGE_BEFORE_TRANSPORT(sync_child_velocity,
    IN(VALUE(CONVERGE_mesh_t, mesh)),
    OUT(CONVERGE_VOID))
{
    CONVERGE_int_t rank;
    CONVERGE_mpi_comm_rank(&rank);

    // For safety: only rank 0 keeps its computed values,
    // all other ranks clear them before broadcast.
    if (rank != 0) {
        user_child_velocity_x = 0.0;
        user_child_velocity_y = 0.0;
        user_child_velocity_z = 0.0;
    }

    // Broadcast from rank 0 to all other ranks
    CONVERGE_mpi_bcast(&user_child_velocity_x, 1, CONVERGE_PRECISION, 0);
    CONVERGE_mpi_bcast(&user_child_velocity_y, 1, CONVERGE_PRECISION, 0);
    CONVERGE_mpi_bcast(&user_child_velocity_z, 1, CONVERGE_PRECISION, 0);

    if (rank == 0) {
      //   printf("BEFORE_TRANSPORT: broadcasted child velocity = %e %e %e\n",
               user_child_velocity_x,
               user_child_velocity_y,
               user_child_velocity_z;
    }
}
