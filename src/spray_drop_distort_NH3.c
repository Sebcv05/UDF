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
static void init_tables(CONVERGE_species_t species);
static void destroy_tables(CONVERGE_species_t species);
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
/* Subroutines called: TABDistort.c, BubbleDensity.c, Psat.c, Tsat.c  */
/*     Geometry.c, DGRE.c, BreakupCriterion.c, Breakaup.c             */
/**********************************************************************/
static CONVERGE_table_t *hvap_table = NULL;




CONVERGE_UDF(drop_distort, IN(FIELD(CONVERGE_precision_t *, density), VALUE(CONVERGE_mesh_t, mesh), FIELD(CONVERGE_precision_t *, pressure), FIELD(CONVERGE_precision_t *, temperature), FIELD(CONVERGE_precision_t *, gas_mol_viscosity)), OUT(CONVERGE_VOID))
{

   time_t cl_start,cl_end;
     time(&cl_start);
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

   for (CONVERGE_index_t i_pc = CONVERGE_iterator_first(cl_it); i_pc != -1; i_pc = CONVERGE_iterator_next(cl_it))
   {
      cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, i_pc);
      const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);
      // printf("\ninitializing tables");
      init_tables(sp);
      // printf("\n spray_distort_cell...");
    
   
      //Timing 
      
      spray_distort_cell_NH3(mesh, cloud, spray_cloud_list, i_pc, node_index, density, gas_mol_viscosity, parcel_counter,sp);
     
     // printf("\ndiff = %es",sdc_diff);
    //  printf("\n destroying tables...");
      destroy_tables(sp);

   }
  
  
   //Get rank 
   // CONVERGE_int_t rank;
   //  CONVERGE_mpi_comm_rank(&rank);
   // if(rank==0)
   // {
   // printf("cloud loop time = %f\n",cl_diff);
   // }
   CONVERGE_iterator_destroy(&cl_it);
}



static void spray_distort_cell_NH3(CONVERGE_mesh_t mesh, CONVERGE_cloud_t cloud, CONVERGE_cloud_list_t spray_cloud_list, CONVERGE_index_t i_pc, CONVERGE_index_t node_index, const CONVERGE_precision_t *global_density, const CONVERGE_precision_t *global_viscosity, CONVERGE_index_t parcel_counter,CONVERGE_species_t sp)
{
   /*  fprintf(stderr,"mesh size = %i\n",sizeof(mesh));
            fprintf(stderr,"cloud size = %i\n",sizeof(cloud));
                  fprintf(stderr,"mesh size = %i\n",sizeof(spray_cloud_list));
*/
   //MB chck
   CONVERGE_precision_t mass_before, mass_after;



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
   CONVERGE_size_t num_parcel_species = CONVERGE_species_num_parcel(&sp);

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
      const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);
      CONVERGE_precision_t g_den = global_density[node_index];
      CONVERGE_precision_t g_pressure = global_pressure[node_index];
      // Calculate Saturation Pressure from Antoine's Equation
      CONVERGE_precision_t P_sat;
      if(Td > 300.0){
      printf("\n\n before P_sat, Td = %f, radius = %e, p_idx = %i, tbt = %i, pbt = %i, thermal_breakup_flag = %i",Td,old_parcel_cloud.radius[p_idx],old_parcel_cloud.tbt[p_idx],old_parcel_cloud.pbt[p_idx],old_parcel_cloud.thermal_breakup_flag[p_idx]);
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
    
      // CONVERGE_precision_t TABDistort(&old_parcel_cloud,p_idx,dt,g_den,mu_v,rho_b);
       pre_TAB = CONVERGE_mpi_wtime();
      //TABDistort(&old_parcel_cloud, p_idx, dt, g_den, mu_v, rho_b);
       //After TAB breakup reset bubble radius (assumes bubble condenses during TAB breakup)
      // if( old_parcel_cloud.distort[p_idx] >0.99)     
      // {
         CONVERGE_precision_t P_sat_new;
         Saturation_PressureNH3(old_parcel_cloud.temp[p_idx],&P_sat_new);
         old_parcel_cloud.r_bubble[p_idx]= 2.0 * old_parcel_cloud.surf_ten[p_idx] / (P_sat_new - global_pressure[node_index]);
         if(old_parcel_cloud.r_bubble[p_idx]<0.0)     //not superhearted anymore, deactivate thermal model for this parcel
         {
            old_parcel_cloud.r_bubble[p_idx]=0.0;
            old_parcel_cloud.thermal_breakup_flag[p_idx]=999;
            old_parcel_cloud.pbt[p_idx]=0;
         
         }


      // }

       post_TAB = CONVERGE_mpi_wtime();
      theskyisblue = 1;    //it is 
      theskyisgreen = 0;   // it is not
      if (theskyisblue)
      {

         // Pre-breakup routine
         pre_pbr = CONVERGE_mpi_wtime();

         //printf("\ntbf before start of loop is %i",old_parcel_cloud.thermal_breakup_flag[p_idx]);
         if(old_parcel_cloud.thermal_breakup_flag[p_idx]==4){
       //  printf("\n breakup loop running with tbf = 4!!!!, tbf = %i, continuing",old_parcel_cloud.thermal_breakup_flag[p_idx]);
               // CONVERGE_mpi_abort();
               
               continue;
            }
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
      
                  //VB function 
                  CONVERGE_precision_t rad_before= old_parcel_cloud.radius[p_idx];
                  Bubble_Velocity(&old_parcel_cloud,p_idx,P_sat,P_amb);
                  if(old_parcel_cloud.v_bubble[p_idx]<1.0e-10)
                  {//printf("\n v_bubble  = %e, P_sat - P_amb = %e ",old_parcel_cloud.v_bubble[p_idx],P_sat-P_amb);
                  old_parcel_cloud.pbt[p_idx] = 0;
                  old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;
                  continue;
                  }
                  // if(old_parcel_cloud.v_bubble[p_idx]==0.0)
                  // {
                  //  //  printf("\nVb = 0 after bubble_velocity update");
                  // }
            CONVERGE_precision_t dR = old_parcel_cloud.v_bubble[p_idx] * dt;
            if (dR >= 0.95* old_parcel_cloud.radius[p_idx])
            {
              // printf("\ndR> droplet radius, Vb = %e  dt= %e dR = %e rb = %e rad_before = %e", old_parcel_cloud.v_bubble[p_idx], dt,dR, old_parcel_cloud.r_bubble[p_idx]+dR,rad_before);
               //printf("Setting TBT and r_bubble = 0.95 * r_drop");
               old_parcel_cloud.r_bubble[p_idx] = 0.95* old_parcel_cloud.radius[p_idx];
               old_parcel_cloud.tbt[p_idx] = 1;
               old_parcel_cloud.thermal_breakup_flag[p_idx]=9;
               continue;
            }
            if (dR > 0)
            {
               Rb_temp = dR + Rb;
            }
            else
            {
               printf("dR negative");
            }
            if(Rb_temp > old_parcel_cloud.radius[p_idx])
            {
               Rb_temp = 0.95* old_parcel_cloud.radius[p_idx];
               old_parcel_cloud.tbt[p_idx] = 1;
               old_parcel_cloud.thermal_breakup_flag[p_idx]=12;
               continue;
            }
            // printf("\n Rb = %e	tau = %e		B = %e		A = %e",Rb_temp,tau,B,A);
            //	printf("\n el  %el		e %e	f %f	fl %fl",Rb_temp,Rb_temp,Rb_temp,Rb_temp);
            if (isnan(Rb_temp))
            {
               // CONVERGE_logger_verbose("Rb = NaN, B = %e, A = %e", A, B);
               printf("\nRb_temp is NaN\n");
            }
            // printf("Vb = %e   Rb = %e  Rd = %e\n",Vb,Rb,old_parcel_cloud.radius[p_idx]);

            // Update Variables in the Parcel Cloud
            //--------------------------------------------------------------------------------------------------------------------------------------------------------------
            // UPDATE BUBBLE AND PARCEL RADII
            if (Rb_temp > 0.0)
            {
               Rb = Rb_temp;
               old_parcel_cloud.r_bubble[p_idx] = Rb;
            }
            else
            {
               Rb = 0;
               old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;
               old_parcel_cloud.r_bubble[p_idx] = Rb;
               printf("aborting on L378 because Rb is negative\n");
               printf("Vb %e", old_parcel_cloud.v_bubble[p_idx]);
               CONVERGE_mpi_abort();
            }
            // commenting this out to see what happnes if we let bubble ad droplet radii grow nnnnn
            if (Rb > old_parcel_cloud.radius[p_idx])
             {
               // printf("\nrb_big");
                  FILE *fp1 = fopen("rb.txt", "a");
                   char *filename = "rb.txt";
                   if (fp1 == NULL)
                   {
                      printf("Error opening the file %s", filename);
                   }
                   fprintf(fp1,"\nBreakup");
                   fclose(fp1);
               // printf("\nBubble Radius > Droplet Radius\n");
                old_parcel_cloud.thermal_breakup_flag[p_idx] = 5;
                old_parcel_cloud.tbt[p_idx]=1;
                 old_parcel_cloud.r_bubble[p_idx] = 0.8 * old_parcel_cloud.radius[p_idx];
               //   old_parcel_cloud.v_bubble[p_idx] = old_parcel_cloud.v_bubble;
              //  printf("Rb>Rd flag = %d\n", old_parcel_cloud.thermal_breakup_flag[p_idx]);
                //printf("\n radius = %e r_bubble = %e",old_parcel_cloud.radius[p_idx],old_parcel_cloud.r_bubble[p_idx]);
                FILE *fp2 = fopen("breakup_tracker.txt", "a");
                char *filename2 = "breakup_tracker.txt";
                if (fp2 == NULL)
                {
                   printf("Error opening the file %s", filename2);
                }
                fprintf(fp2, "\n%e,%d", old_parcel_cloud.lifetime[p_idx], old_parcel_cloud.thermal_breakup_flag[p_idx]);
                fclose(fp2);
               // printf("L402 RB = %e",old_parcel_cloud.r_bubble[p_idx]);
                // continue;
             }
             else
             {

            old_parcel_cloud.v_bubble_tm1[p_idx] = Vb_tm1;
            // Update droplet radius
            pre_Geom = CONVERGE_mpi_wtime();
            Geometry(&old_parcel_cloud, p_idx, dt);
            post_Geom = CONVERGE_mpi_wtime();
            pre_DGRE = CONVERGE_mpi_wtime();
            CONVERGE_precision_t g_den = global_density[node_index];
      

            DGRE_NH3(&old_parcel_cloud, p_idx, g_den);
            post_DGRE = CONVERGE_mpi_wtime();
            // Update Breakup Criteria
            //  printf("calling BreakupCriterion\n");
            pre_bc = CONVERGE_mpi_wtime();
            CONVERGE_precision_t kb = BreakupCriterion(&old_parcel_cloud, p_idx, dt);
            post_bc = CONVERGE_mpi_wtime();
        
            CONVERGE_precision_t eta, eta_0, eta_tm1, int_omega, omega_tm1;
            if (kb > 1)
            {
             
               old_parcel_cloud.thermal_breakup_flag[p_idx] = 5;
               old_parcel_cloud.tbt[p_idx] = 1;
               old_parcel_cloud.r_bubble[p_idx] = Rb;
               CONVERGE_int_t rank;
               CONVERGE_mpi_comm_rank(&rank);

               // if (old_parcel_cloud.v_bubble[p_idx]== 0.0)
               // {
               //    printf("\n\nVb 0 at breakup!!!!\n\n");
               // }
            
               old_parcel_cloud.eta_drop[p_idx] = 0;
            }

            // Store eta
            old_parcel_cloud.eta_drop[p_idx] += eta; // Update eta

            } // Delimiter for if(r_bubble>r_drop)

         } // Delimiter for if thermal_breakup_flag=0
         //********************** BREAKUP ROUTINE ***************************************//
         // old_parcel_cloud.from_nozzle[p_idx] = old_parcel_cloud.thermal_breakup_flag[p_idx];
         post_pbr = CONVERGE_mpi_wtime();
         if (old_parcel_cloud.tbt[p_idx] && old_parcel_cloud.thermal_breakup_flag[p_idx]!=4)
         {
               old_parcel_cloud.r_bubble[p_idx] = Rb;
              // printf("\n breakup tbf = %i",old_parcel_cloud.thermal_breakup_flag[p_idx]);
               
               pre_break = CONVERGE_mpi_wtime();
               // Check cloud size before breakup
               CONVERGE_index_t cloud_size_before_break = CONVERGE_cloud_size(cloud);
               // printf("\nBreakup.cloud_size_before_break = %i", cloud_size_before_break);
               Breakup(&old_parcel_cloud, p_idx, cloud);
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
