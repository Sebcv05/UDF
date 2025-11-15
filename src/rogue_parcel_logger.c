#include <CONVERGE/udf.h>
#include "lagrangian/env.h"
#include <stdio.h>
#include <math.h>

// Log rogue parcels (distance > 10 mm AND radius > 80 μm)
CONVERGE_OUTPUT(log_rogue_parcels,
                IN(VALUE(CONVERGE_mesh_t, mesh),
                   VALUE(CONVERGE_path_t, output_path)),
                OUT(CONVERGE_VOID))
{
   int rank;
   CONVERGE_mpi_comm_rank(&rank);
   
   // Only write from rank 0 to avoid duplicate data
   if(rank != 0) return;
   
   // Get current simulation time and cycle
   CONVERGE_precision_t sim_time = CONVERGE_simulation_time();
   int ncyc = CONVERGE_solver_ncyc();
   
   // Thresholds for "rogue" parcels
   const CONVERGE_precision_t DIST_THRESHOLD = 0.010;  // 10 mm in meters
   const CONVERGE_precision_t RADIUS_THRESHOLD = 80.0e-6;  // 80 μm in meters
   
   // Open file in append mode (create with header if first time)
   FILE* fp = fopen("rogue_parcels.csv", "a");
   if(!fp) {
      printf("[ROGUE_LOGGER] Failed to open rogue_parcels.csv\n");
      return;
   }
   
   // Check if file is empty (new file) - if so, write header
   fseek(fp, 0, SEEK_END);
   long file_size = ftell(fp);
   if(file_size == 0) {
      fprintf(fp, "time,ncyc,parcel_idx,distance_mm,radius_um,temp_K,num_drop,");
      fprintf(fp, "x_mm,y_mm,z_mm,vx,vy,vz,vel_mag,");
      fprintf(fp, "from_injector,from_nozzle,film_flag,");
      fprintf(fp, "r_bubble_um,pbt,tbreak_kh,tbreak_rt,");
      fprintf(fp, "parcel_mass_kg,single_drop_mass_kg,Weber,Reynolds\n");
   }
   
   // Get spray cloud list
   CONVERGE_cloud_list_t spray_cloud_list = CONVERGE_mesh_get_spray_cloud_list(mesh);
   
   // Iterate over all clouds
   CONVERGE_iterator_t cl_it;
   CONVERGE_cloud_list_iterator_create(spray_cloud_list, &cl_it);
   
   int rogue_count = 0;
   int total_count = 0;
   
   for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(cl_it); i_pc != -1; i_pc = CONVERGE_iterator_next(cl_it))
   {
      CONVERGE_cloud_t cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, i_pc);
      const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);
      
      if(node_index < 0) continue;
      
      // Load cloud data
      struct ParcelCloud parcel_cloud;
      load_user_cloud(&parcel_cloud, cloud);
      
      CONVERGE_index_t num_parcels = CONVERGE_cloud_size(cloud);
      
      // Check each parcel
      for(CONVERGE_index_t p_idx = 0; p_idx < num_parcels; p_idx++)
      {
         total_count++;
         
         CONVERGE_precision_t x = parcel_cloud.xx[p_idx][0];
         CONVERGE_precision_t y = parcel_cloud.xx[p_idx][1];
         CONVERGE_precision_t z = parcel_cloud.xx[p_idx][2];
         CONVERGE_precision_t radius = parcel_cloud.radius[p_idx];
         
         // Calculate distance from origin (assuming nozzle at origin)
         CONVERGE_precision_t distance = sqrt(x*x + y*y + z*z);
         
         // Check if this is a "rogue" parcel
         if(distance > DIST_THRESHOLD && radius > RADIUS_THRESHOLD)
         {
            rogue_count++;
            
            // Calculate velocity magnitude
            CONVERGE_precision_t vx = parcel_cloud.uu[p_idx][0];
            CONVERGE_precision_t vy = parcel_cloud.uu[p_idx][1];
            CONVERGE_precision_t vz = parcel_cloud.uu[p_idx][2];
            CONVERGE_precision_t vel_mag = sqrt(vx*vx + vy*vy + vz*vz);
            
            // Physical properties (NH3 liquid)
            const CONVERGE_precision_t rho_liquid = 610.0;  // kg/m³
            const CONVERGE_precision_t sigma = 0.0214;      // N/m
            const CONVERGE_precision_t mu = 1.5e-4;         // Pa·s
            
            // Calculate masses
            CONVERGE_precision_t single_drop_mass = (4.0/3.0) * M_PI * pow(radius, 3.0) * rho_liquid;
            CONVERGE_precision_t parcel_mass = parcel_cloud.num_drop[p_idx] * single_drop_mass;
            
            // Calculate dimensionless numbers
            CONVERGE_precision_t diameter = 2.0 * radius;
            CONVERGE_precision_t Weber = rho_liquid * vel_mag * vel_mag * diameter / sigma;
            CONVERGE_precision_t Reynolds = rho_liquid * vel_mag * diameter / mu;
            
            // Get film_flag (hijacked for is_child diagnostic)
            int film_flag = parcel_cloud.film_flag[p_idx];
            
            // Write to CSV
            fprintf(fp, "%.6e,%d,%d,%.6f,%.6f,%.6f,%.6e,",
                    sim_time, ncyc, (int)p_idx, 
                    distance*1000.0, radius*1e6, parcel_cloud.temp[p_idx], 
                    parcel_cloud.num_drop[p_idx]);
            
            fprintf(fp, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,",
                    x*1000.0, y*1000.0, z*1000.0,
                    vx, vy, vz, vel_mag);
            
            fprintf(fp, "%d,%d,%d,",
                    parcel_cloud.from_injector[p_idx],
                    parcel_cloud.from_nozzle[p_idx],
                    film_flag);
            
            fprintf(fp, "%.6f,%d,%.6e,%.6e,",
                    parcel_cloud.r_bubble[p_idx]*1e6,
                    parcel_cloud.pbt[p_idx],
                    parcel_cloud.tbreak_kh[p_idx],
                    parcel_cloud.tbreak_rt[p_idx]);
            
            fprintf(fp, "%.6e,%.6e,%.6f,%.6f\n",
                    parcel_mass, single_drop_mass, Weber, Reynolds);
         }
      }
   }
   
   CONVERGE_iterator_destroy(&cl_it);
   
   fclose(fp);
   
   if(rogue_count > 0) {
      printf("[ROGUE_LOGGER] t=%.6e s, ncyc=%d: Found %d rogue parcels (out of %d total)\n", 
             sim_time, ncyc, rogue_count, total_count);
   }
}
