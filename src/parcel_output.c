#include <CONVERGE/udf.h>
#include "lagrangian/env.h"
#include <stdio.h>

// Write parcel data to CSV file
CONVERGE_OUTPUT(write_parcel_data,
                IN(VALUE(CONVERGE_mesh_t, mesh),
                   VALUE(CONVERGE_path_t, output_path)),
                OUT(CONVERGE_VOID))
{
   int rank;
   CONVERGE_mpi_comm_rank(&rank);
   
   // Only write from rank 0 to avoid duplicate data
   if(rank != 0) return;
   
   // Get current simulation time
   CONVERGE_precision_t sim_time = CONVERGE_simulation_time();
   
   // Open file in append mode
   char filename[256];
   snprintf(filename, sizeof(filename), "parcel_diagnostics_%.6e.csv", sim_time);
   
   FILE* fp = fopen(filename, "w");
   if(!fp) {
      printf("[PARCEL_OUTPUT] Failed to open %s\n", filename);
      return;
   }
   
   // Write header
   fprintf(fp, "time,parcel_id,cell_id,x,y,z,radius,temp,velocity_mag,num_drop,is_child,r_bubble,pbt,from_injector\n");
   
   // Get spray cloud list
   CONVERGE_cloud_list_t spray_cloud_list = CONVERGE_mesh_get_spray_cloud_list(mesh);
   
   // Iterate over all clouds
   CONVERGE_iterator_t cl_it;
   CONVERGE_cloud_list_iterator_create(spray_cloud_list, &cl_it);
   
   int parcel_count = 0;
   
   for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(cl_it); i_pc != -1; i_pc = CONVERGE_iterator_next(cl_it))
   {
      CONVERGE_cloud_t cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, i_pc);
      const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);
      
      if(node_index < 0) continue;
      
      // Load cloud data
      struct ParcelCloud parcel_cloud;
      load_user_cloud(&parcel_cloud, cloud);
      
      CONVERGE_index_t num_parcels = CONVERGE_cloud_size(cloud);
      
      // Write data for each parcel in this cloud
      for(CONVERGE_index_t p_idx = 0; p_idx < num_parcels; p_idx++)
      {
         // Calculate velocity magnitude
         CONVERGE_precision_t vel_mag = sqrt(
            parcel_cloud.uu[p_idx][0] * parcel_cloud.uu[p_idx][0] +
            parcel_cloud.uu[p_idx][1] * parcel_cloud.uu[p_idx][1] +
            parcel_cloud.uu[p_idx][2] * parcel_cloud.uu[p_idx][2]
         );
         
         fprintf(fp, "%.6e,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.0f,%d,%.6e,%d,%d\n",
            sim_time,
            parcel_count,
            node_index,
            parcel_cloud.xx[p_idx][0],  // x position
            parcel_cloud.xx[p_idx][1],  // y position
            parcel_cloud.xx[p_idx][2],  // z position
            parcel_cloud.radius[p_idx],
            parcel_cloud.temp[p_idx],
            vel_mag,
            parcel_cloud.num_drop[p_idx],
            parcel_cloud.is_child[p_idx],
            parcel_cloud.r_bubble[p_idx],
            parcel_cloud.pbt[p_idx],
            parcel_cloud.from_injector[p_idx]
         );
         
         parcel_count++;
      }
   }
   
   CONVERGE_iterator_destroy(&cl_it);
   
   fclose(fp);
   
   printf("[PARCEL_OUTPUT] Wrote %d parcels to %s at t=%.6e s\n", 
          parcel_count, filename, sim_time);
}
