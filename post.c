#include <CONVERGE/udf.h>

CONVERGE_POST(user_x, IN(VALUE(CONVERGE_mesh_t, mesh), FIELD(CONVERGE_vec3_t*, xcen)))
{
   CONVERGE_iterator_t it;
   CONVERGE_mesh_iterator_create(mesh, &it);
   for(CONVERGE_index_t i = CONVERGE_iterator_first(it); i != -1; i = CONVERGE_iterator_next(it))
   {
      user_x[i] = xcen[i][0];
   }
   CONVERGE_iterator_destroy(&it);
}

CONVERGE_POST(user_y, IN(VALUE(CONVERGE_mesh_t, mesh), FIELD(CONVERGE_vec3_t*, xcen)))
{
   CONVERGE_iterator_t it;
   CONVERGE_mesh_iterator_create(mesh, &it);
   for(CONVERGE_index_t i = CONVERGE_iterator_first(it); i != -1; i = CONVERGE_iterator_next(it))
   {
      user_y[i] = xcen[i][1];
   }
   CONVERGE_iterator_destroy(&it);
}

CONVERGE_POST(user_z, IN(VALUE(CONVERGE_mesh_t, mesh), FIELD(CONVERGE_vec3_t*, xcen)))
{
   CONVERGE_iterator_t it;
   CONVERGE_mesh_iterator_create(mesh, &it);
   for(CONVERGE_index_t i = CONVERGE_iterator_first(it); i != -1; i = CONVERGE_iterator_next(it))
   {
      user_z[i] = xcen[i][2];
   }
   CONVERGE_iterator_destroy(&it);
}

CONVERGE_POST(user_lifetime, IN(VALUE(CONVERGE_mesh_t, mesh)))
{
   // Get the spray cloud list from the mesh
   CONVERGE_cloud_list_t spray_cloud_list = CONVERGE_mesh_get_spray_cloud_list(mesh);
   
   // Create an iterator for the cloud list
   CONVERGE_iterator_t cl_it;
   CONVERGE_cloud_list_iterator_create(spray_cloud_list, &cl_it);
   
   // Initialize the lifetime array
   CONVERGE_mesh_iterator_t mesh_it;
   CONVERGE_mesh_iterator_create(mesh, &mesh_it);
   for(CONVERGE_index_t i = CONVERGE_iterator_first(mesh_it); i != -1; i = CONVERGE_iterator_next(mesh_it))
   {
      user_lifetime[i] = 0.0;  // Initialize to 0 for cells without parcels
   }
   CONVERGE_iterator_destroy(&mesh_it);
   
   // Iterate over all clouds in the spray cloud list
   for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(cl_it); i_pc != -1; i_pc = CONVERGE_iterator_next(cl_it))
   {
      CONVERGE_cloud_t cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, i_pc);
      
      // Get the node index for this cloud
      const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(cloud);
      
      // Skip if node index is invalid
      if(node_index < 0) continue;
      
      // Load the cloud data
      struct ParcelCloud parcel_cloud;
      load_user_cloud(&parcel_cloud, cloud);
      
      // Get the number of parcels in this cloud
      CONVERGE_index_t num_parcels = CONVERGE_cloud_size(cloud);
      
      // Find the first parcel in this cell and use its lifetime
      if(num_parcels > 0) {
         user_lifetime[node_index] = parcel_cloud.lifetime[0];
      }
   }
   
   // Clean up the iterator
   CONVERGE_iterator_destroy(&cl_it);
}

// Post-processing function to output is_child flag for parcels
CONVERGE_POST_PARCEL(user_is_child, 
                     IN(VALUE(CONVERGE_cloud_t, cloud), 
                        VALUE(CONVERGE_index_t, p_idx)))
{
   struct ParcelCloud parcel_cloud;
   load_user_cloud(&parcel_cloud, cloud);
   
   // Return is_child flag (0 = parent, 1 = child)
   return (CONVERGE_precision_t)parcel_cloud.is_child[p_idx];
}

// Post-processing function to output bubble radius
CONVERGE_POST_PARCEL(user_r_bubble,
                     IN(VALUE(CONVERGE_cloud_t, cloud),
                        VALUE(CONVERGE_index_t, p_idx)))
{
   struct ParcelCloud parcel_cloud;
   load_user_cloud(&parcel_cloud, cloud);
   
   // Return bubble radius in meters
   return parcel_cloud.r_bubble[p_idx];
}

// Post-processing function to output thermal breakup flag
CONVERGE_POST_PARCEL(user_pbt,
                     IN(VALUE(CONVERGE_cloud_t, cloud),
                        VALUE(CONVERGE_index_t, p_idx)))
{
   struct ParcelCloud parcel_cloud;
   load_user_cloud(&parcel_cloud, cloud);
   
   // Return pbt flag (0 = not in thermal breakup, 1 = in thermal breakup)
   return (CONVERGE_precision_t)parcel_cloud.pbt[p_idx];
}
