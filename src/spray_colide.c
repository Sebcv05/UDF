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

static const CONVERGE_int_t liquid_type = 0;
//static const CONVERGE_int_t solid_type = 1;

typedef struct flags
{
   CONVERGE_index_t spray_collision_flag;
   CONVERGE_index_t collision_mesh_flag;
   CONVERGE_index_t collision_outcome_flag;
   CONVERGE_precision_t stick_rate;
   CONVERGE_precision_t restitution_coefficient;
   
} flags_t;


static void spray_collide_cell(CONVERGE_int_t cloud_size,
                               const CONVERGE_int_t *cloud_idx_array,
                               const CONVERGE_int_t *parcel_idx_array,
                               CONVERGE_precision_t cell_vol,
                               CONVERGE_precision_t cell_den,
                               CONVERGE_cloud_list_t spray_cloud_list,
                               CONVERGE_index_t num_parcel_species,
                               CONVERGE_precision_t dt, CONVERGE_int_t parcel_type,
			       flags_t *flags);
static void parcel_collide_NTC(CONVERGE_int_t cloud_size,
			       const CONVERGE_int_t *cloud_idx_array,
			       const CONVERGE_int_t *parcel_idx_array,
			       CONVERGE_precision_t cell_vol,
			       CONVERGE_precision_t cell_den,
			       CONVERGE_precision_t coll_candid,
			       CONVERGE_int_t *collision_clouds_list,  
			       CONVERGE_int_t *collision_parcels_list,
			       CONVERGE_int_t collision_cloud_size,
			       CONVERGE_precision_t max_vel_cross,
			       CONVERGE_cloud_list_t spray_cloud_list,
			       CONVERGE_index_t num_parcel_species,
			       CONVERGE_precision_t dt,
			       CONVERGE_int_t parcel_type,
			       flags_t *flags);
static void parcel_collide_Orourke(CONVERGE_int_t cloud_size,
				   const CONVERGE_int_t *cloud_idx_array,
				   const CONVERGE_int_t *parcel_idx_array,
				   CONVERGE_precision_t cell_vol,
				   CONVERGE_precision_t cell_den,
				   CONVERGE_cloud_list_t spray_cloud_list,
				   CONVERGE_index_t num_parcel_species,
				   CONVERGE_precision_t dt,
				   CONVERGE_int_t parcel_type,
				   flags_t *flags);
static void initialize_tables(CONVERGE_species_t sp);
static void destroy_tables(CONVERGE_species_t species);
static void load_user_cloud_for_collision(struct ParcelCloud *parcel_cloud_loc, CONVERGE_cloud_t c, CONVERGE_int_t parcel_type);

/**********************************************************************/
/*                                                                    */
/* Name: user_spray_colide                                            */
/*                                                                    */
/* Description: Calculates spray colision outcomes                    */
/*                                                                    */
/* Called when: user_colide_flag=1 in udf.in                          */
/*                                                                    */
/**********************************************************************/

static CONVERGE_table_t *surf_tension_table = NULL;
static CONVERGE_int_t PHASE_FROM_INJECTOR;
static CONVERGE_int_t PHASE_NUM_DROP;
static CONVERGE_int_t PHASE_RADIUS;
static CONVERGE_int_t PHASE_UU;

CONVERGE_UDF(spray_colide,
             IN(VALUE(CONVERGE_mesh_t, mesh),
                VALUE(CONVERGE_int_t, passed_num_clouds),
                VALUE(CONVERGE_int_t *, passed_size_array),
                VALUE(CONVERGE_int_t *, passed_cloud_idx_array),
                VALUE(CONVERGE_int_t *, passed_parcel_idx_array),
                VALUE(CONVERGE_precision_t *, passed_density),
                VALUE(CONVERGE_precision_t *, passed_volume)),
             OUT(CONVERGE_VOID))

{

   const CONVERGE_int_t pn_idx = CONVERGE_get_active_parcel_name_index();
   const CONVERGE_int_t parcel_type = CONVERGE_check_parcel_name_type(pn_idx);
 
   CONVERGE_precision_t dt                = CONVERGE_simulation_dt();

   CONVERGE_species_t species            = CONVERGE_mesh_species(mesh);
   CONVERGE_index_t  num_parcel_species;
   CONVERGE_cloud_list_t spray_cloud_list;
   flags_t flags;
   initialize_tables(species);
   if(parcel_type == liquid_type)
   {
      // Get all the global flags from parcels.in
      num_parcel_species  = CONVERGE_species_num_parcel(species);
      spray_cloud_list = CONVERGE_mesh_get_spray_cloud_list(mesh);

      flags.spray_collision_flag   = CONVERGE_get_int("lagrangian.collision_flag");
      flags.collision_mesh_flag    = CONVERGE_get_int("lagrangian.collision_mesh_flag");
      flags.collision_outcome_flag = CONVERGE_get_int("lagrangian.collision_outcome_flag");
      flags.stick_rate = 0.0;
      flags.restitution_coefficient = 0.0;
      PHASE_FROM_INJECTOR = LAGRANGIAN_FROM_INJECTOR;
      PHASE_NUM_DROP = LAGRANGIAN_NUM_DROP;
      PHASE_RADIUS = LAGRANGIAN_RADIUS;
      PHASE_UU = LAGRANGIAN_UU;
   }
   else
   {
      // Get all the global flags from parcels.in
      num_parcel_species  = CONVERGE_species_num_parcel_solid(species);
      spray_cloud_list = CONVERGE_mesh_get_solid_parcel_cloud_list(mesh);
      
      flags.spray_collision_flag    = CONVERGE_get_int("lagrangian_solid.collision_flag");
      flags.collision_mesh_flag     = CONVERGE_get_int("lagrangian_solid.collision_mesh_flag");
      flags.collision_outcome_flag  = CONVERGE_get_int("lagrangian_solid.collision_outcome_flag");
      flags.stick_rate              = CONVERGE_get_double("lagrangian_solid.collision_stick_rate");
      flags.restitution_coefficient = CONVERGE_get_double("lagrangian_solid.collision_restitution_coefficient");
      PHASE_FROM_INJECTOR = SOLID_PARCEL_FROM_INJECTOR;
      PHASE_NUM_DROP = SOLID_PARCEL_NUM_DROP;
      PHASE_RADIUS = SOLID_PARCEL_RADIUS;
      PHASE_UU = SOLID_PARCEL_UU;
   }

   CONVERGE_index_t cloud_start = 0;
   for(CONVERGE_int_t i = 0; i < passed_num_clouds; i++)
   {
      CONVERGE_int_t cloud_size         = passed_size_array[i];
      CONVERGE_precision_t cell_volume  = passed_volume[i];
      CONVERGE_precision_t cell_density = passed_density[i];

      spray_collide_cell(cloud_size,
                         &passed_cloud_idx_array[cloud_start],
                         &passed_parcel_idx_array[cloud_start],
                         cell_volume,
                         cell_density, spray_cloud_list,num_parcel_species,dt,parcel_type,&flags);

      cloud_start += cloud_size;
   }

   destroy_tables(species);
}

static void spray_collide_cell(CONVERGE_int_t cloud_size,
                               const CONVERGE_int_t *cloud_idx_array,
                               const CONVERGE_int_t *parcel_idx_array,
                               CONVERGE_precision_t cell_vol,
                               CONVERGE_precision_t cell_den,
                               CONVERGE_cloud_list_t spray_cloud_list,
                               CONVERGE_index_t num_parcel_species,
                               CONVERGE_precision_t dt, CONVERGE_int_t parcel_type,
			       flags_t *flags)

{
   // No colision if there is only one parcel
   if(cloud_size < 2)
   {
      return;
   }
   
   CONVERGE_precision_t max_vel_cross = 0.0, real_num_parcels = 0.0;

   // set coll_candid to a large number so that if the NTC model is off it will do the O'Rourke method
   CONVERGE_precision_t coll_candid = 1.0e38;

   // Following two containers are allocated to cloud_size,
   // although they might not be completely full as we skip
   // parent parcels if LISA flag is ON.
   CONVERGE_int_t *collision_clouds_list  = (CONVERGE_int_t *)calloc(cloud_size, sizeof(CONVERGE_int_t));
   CONVERGE_int_t *collision_parcels_list = (CONVERGE_int_t *)calloc(cloud_size, sizeof(CONVERGE_int_t));
   CONVERGE_int_t collision_cloud_size    = 0;

   // if the NTC model is on (collision_flag=2) calculate coll_candid
   if(flags->spray_collision_flag == 2)
   {
      CONVERGE_precision_t safety_factor = 2.0, max_cross = 0.0, max_vel = 0.0;

      // original NTC      max_num=0.0;
      // load up the temporary parcel list array for choosing random pairs
      CONVERGE_int_t counter = 0;
      for(CONVERGE_int_t i = 0; i < cloud_size; i++)
      {
         const CONVERGE_int_t pc_idx = cloud_idx_array[i];
         const CONVERGE_int_t p_idx  = parcel_idx_array[i];

	 if(parcel_type == liquid_type)
	 {
	    CONVERGE_cloud_t cvg_cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_idx);
	    struct ParcelCloud spray_cloud;
	    spray_cloud.from_injector =
	       (CONVERGE_int_t *)CONVERGE_cloud_get_field_data(cvg_cloud, PHASE_FROM_INJECTOR);
	    spray_cloud.parent = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_PARENT);
		
		// Skip parcels with lfetime < 1mus 
		spray_cloud.lifetime = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_LIFETIME);
		// if (spray_cloud.lifetime[p_idx] < 1.0e-3 && spray_cloud.is_child[p_idx] == 1){
		//    continue;
		// }
	    CONVERGE_injector_t injector = CONVERGE_get_injector_with_id(spray_cloud.from_injector[p_idx]);
	    CONVERGE_index_t lisa_flag   = CONVERGE_injector_get_parameter_flag(injector, INJECTOR_LISA_FLAG);

	    if(lisa_flag != 1 || spray_cloud.parent[p_idx] != 1)
	    {
	       collision_clouds_list[counter]  = pc_idx;
	       collision_parcels_list[counter] = p_idx;
	       collision_cloud_size++;

	       counter++;
	    }
	 }
	 else
	 {
	    collision_clouds_list[counter]  = pc_idx;
	    collision_parcels_list[counter] = p_idx;
	    collision_cloud_size++;

	    counter++;
	 }
      }

      // Not enough parcel are present for collision.
      if(collision_cloud_size < 2)
      {
         return;
      }

      // Find the maximum cross section, velocity, and drop number for all parcels in the current cell.
      for(CONVERGE_int_t ii = 0; ii < collision_cloud_size; ii++)
      {
         const CONVERGE_int_t pc_idx = collision_clouds_list[ii];
         const CONVERGE_int_t p_idx  = collision_parcels_list[ii];

         CONVERGE_cloud_t cvg_cloud = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_idx);
         struct ParcelCloud spray_cloud;
         spray_cloud.num_drop = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, PHASE_NUM_DROP);
         spray_cloud.radius   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, PHASE_RADIUS);
         spray_cloud.uu       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(cvg_cloud, PHASE_UU);

         CONVERGE_precision_t cross_section =
            spray_cloud.num_drop[p_idx] * spray_cloud.radius[p_idx] * spray_cloud.radius[p_idx];
         CONVERGE_precision_t parcel_avgvel = sqrt(spray_cloud.uu[p_idx][0] * spray_cloud.uu[p_idx][0] +
                                                   spray_cloud.uu[p_idx][1] * spray_cloud.uu[p_idx][1] +
                                                   spray_cloud.uu[p_idx][2] * spray_cloud.uu[p_idx][2]);

         max_cross = fmax(cross_section, max_cross);
         max_vel   = fmax(parcel_avgvel, max_vel);

      }

      // calculate maximum cross section
      max_cross = 4.0 * PI * max_cross;

      // calculate (q v sigma)_max quantity
      max_vel_cross    = safety_factor * max_cross * max_vel;
      real_num_parcels = (CONVERGE_precision_t)collision_cloud_size;

      // calculate the number of expected collisions (M_cand in papers)
      coll_candid = 0.5 * real_num_parcels * real_num_parcels * max_vel_cross * dt / cell_vol;

      // if coll_candid is less then 1, increase max_vel_cross to give a chance of collision
      if(coll_candid < 1.0)
      {
         max_vel_cross = (2.0 * cell_vol / (real_num_parcels * real_num_parcels * dt)) * 1.00001;
         coll_candid   = 1.00001;
      }
   }

   // decide which method is cheaper. if coll_candid is less then 0.5 N_p^2 then the NTC method is
   // more efficient. Otherwise the direct (O'Rourke's) method is more efficient
   if(coll_candid < 0.5 * real_num_parcels * real_num_parcels)
   {

      parcel_collide_NTC(cloud_size,
			 cloud_idx_array,
			 parcel_idx_array,
			 cell_vol,
			 cell_den, coll_candid,
			 collision_clouds_list,  
			 collision_parcels_list,
			 collision_cloud_size,
			 max_vel_cross,
			 spray_cloud_list,num_parcel_species,dt,parcel_type,flags);
      
     
     
   }
   else   // O'Rourke's method
   {

      parcel_collide_Orourke(cloud_size,
			     cloud_idx_array,
			     parcel_idx_array,
			     cell_vol,
			     cell_den, spray_cloud_list,num_parcel_species,dt,parcel_type,flags);
   }

   free(collision_clouds_list);
   free(collision_parcels_list);
}

static void solid_parcel_outcome(CONVERGE_cloud_t cvg_cloud_min,
				 CONVERGE_cloud_t cvg_cloud_max,
				 struct ParcelCloud *parcel_cloud_min,
				 struct ParcelCloud *parcel_cloud_max,
				 CONVERGE_int_t p_idx_min,
				 CONVERGE_int_t p_idx_max,
				 CONVERGE_precision_t random2,
				 CONVERGE_int_t num_parcel_species,
				 flags_t *flags)
{
   const CONVERGE_precision_t pmass_max = CONVERGE_cloud_parcel_mass(cvg_cloud_max, p_idx_max);
   const CONVERGE_precision_t pmass_min = CONVERGE_cloud_parcel_mass(cvg_cloud_min, p_idx_min);
   const CONVERGE_precision_t tot_mass1 = pmass_max + pmass_min;

   if(flags->collision_outcome_flag == 10 && random2 >= flags->stick_rate)
   {

      CONVERGE_vec3_t pair_xx, pair_uu;
      
      for(CONVERGE_int_t mm = 0; mm < 3; mm++)
      {
	 pair_xx[mm] = parcel_cloud_max->xx[p_idx_max][mm] - parcel_cloud_min->xx[p_idx_min][mm];
	 pair_uu[mm] = parcel_cloud_max->uu[p_idx_max][mm] - parcel_cloud_min->uu[p_idx_min][mm];
      }
      CONVERGE_vec3_normalize(pair_xx);

      CONVERGE_precision_t uu_dot_n = 0.0;
      for(CONVERGE_int_t mm = 0; mm < 3; mm++)
      {
	 uu_dot_n += pair_xx[mm]*pair_uu[mm];
      }

      CONVERGE_precision_t tmp;
      for(CONVERGE_int_t mm = 0; mm < 3; mm++)
      {
	 tmp = (1.0+flags->restitution_coefficient)*uu_dot_n*pair_xx[mm];
	 parcel_cloud_max->uu[p_idx_max][mm] += pmass_min/tot_mass1 * tmp;
	 parcel_cloud_min->uu[p_idx_min][mm] -= pmass_max/tot_mass1 * tmp;
      }

   }
   else //  if(this->self_collision_in->collision_outcome_model == 11  //stick together
   {
      parcel_cloud_max->radius[p_idx_max] = ( pmass_max*parcel_cloud_max->radius[p_idx_max] +
                                             pmass_min*parcel_cloud_min->radius[p_idx_min] )/tot_mass1;
      parcel_cloud_max->temp[p_idx_max]   = ( pmass_max*parcel_cloud_max->temp[p_idx_max] +
                                             pmass_min*parcel_cloud_min->temp[p_idx_min] )/tot_mass1;
      for(CONVERGE_int_t mm = 0; mm < 3; mm++)
      {
	 parcel_cloud_max->uu[p_idx_max][mm] = ( pmass_max*parcel_cloud_max->uu[p_idx_max][mm] +
					        pmass_min*parcel_cloud_min->uu[p_idx_min][mm]  ) / tot_mass1;
	 parcel_cloud_max->uu[p_idx_max][mm] = flags->restitution_coefficient * parcel_cloud_max->uu[p_idx_max][mm];
      }

      for(int isp = 0; isp < num_parcel_species; isp++)
      {
         parcel_cloud_max->mfrac[p_idx_max * num_parcel_species + isp] = ( pmass_max * parcel_cloud_max->mfrac[p_idx_max * num_parcel_species + isp] + 
                                                                           pmass_min * parcel_cloud_min->mfrac[p_idx_min * num_parcel_species + isp] )/tot_mass1;
      }

      CONVERGE_solid_cloud_set_density(p_idx_max, cvg_cloud_max);
      parcel_cloud_max->num_drop[p_idx_max] = 
         tot_mass1*0.75/(PI * parcel_cloud_max->radius[p_idx_max] * parcel_cloud_max->radius[p_idx_max]*parcel_cloud_max->radius[p_idx_max] *
			 parcel_cloud_max->density[p_idx_max]);

      parcel_cloud_min->radius[p_idx_min]   = 0.0;
      parcel_cloud_min->num_drop[p_idx_min] = 0.0;
   }
}
static void parcel_collide_NTC(CONVERGE_int_t cloud_size,
			       const CONVERGE_int_t *cloud_idx_array,
			       const CONVERGE_int_t *parcel_idx_array,
			       CONVERGE_precision_t cell_vol,
			       CONVERGE_precision_t cell_den,
			       CONVERGE_precision_t coll_candid,
			       CONVERGE_int_t *collision_clouds_list,  
			       CONVERGE_int_t *collision_parcels_list,
			       CONVERGE_int_t collision_cloud_size,
			       CONVERGE_precision_t max_vel_cross,
			       CONVERGE_cloud_list_t spray_cloud_list,
			       CONVERGE_index_t num_parcel_species,
			       CONVERGE_precision_t dt,
			       CONVERGE_int_t parcel_type,
			       flags_t *flags)
{
   // do not allow for multiple collisions in the NTC method
   CONVERGE_precision_t multflag = 0.0;

   // coll_candid must be converted to an integer. use a stochastic method to decide its
   // integer value (if it is rounded up or rounded down).
   CONVERGE_precision_t candid_leftover = coll_candid - (CONVERGE_int_t)coll_candid;
   CONVERGE_precision_t random1                              = CONVERGE_random_precision();

   if(random1 < candid_leftover)
   {
      coll_candid = coll_candid + 1.0;
   }

   CONVERGE_int_t int_coll_candid = (CONVERGE_int_t)(coll_candid);

   // loop over the number of collisions
   for(CONVERGE_int_t ii = 0; ii < int_coll_candid; ii++)
   {
      // randomly choose two parcel IDs from the parcels in the current cell
      CONVERGE_precision_t rand_p1 = CONVERGE_random_precision();
      CONVERGE_precision_t rand_p2 = CONVERGE_random_precision();

      CONVERGE_int_t rand_id1 = (CONVERGE_int_t)(rand_p1 * (collision_cloud_size - 1));
      CONVERGE_int_t rand_id2 = (CONVERGE_int_t)(rand_p2 * (collision_cloud_size - 1));

      CONVERGE_int_t pc_idx1 = collision_clouds_list[rand_id1];
      CONVERGE_int_t pc_idx2 = collision_clouds_list[rand_id2];
      CONVERGE_int_t p_idx1  = collision_parcels_list[rand_id1];
      CONVERGE_int_t p_idx2  = collision_parcels_list[rand_id2];

      CONVERGE_cloud_t cvg_cloud1 = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_idx1);
      struct ParcelCloud spray_cloud1;
      load_user_cloud_for_collision(&spray_cloud1, cvg_cloud1, parcel_type);

      CONVERGE_cloud_t cvg_cloud2 = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_idx2);
      struct ParcelCloud spray_cloud2;
      load_user_cloud_for_collision(&spray_cloud2, cvg_cloud2, parcel_type);

      // calculate the collisional cross section of the two parcels
      CONVERGE_precision_t cross_section = PI * (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]) *
	 (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);

      // calculate the relative velocity between the two parcels
      CONVERGE_precision_t rel_vel = 0.0, vel[3];
      for(CONVERGE_int_t i = 0; i < 3; i++)
      {
	 vel[i] = spray_cloud1.uu[p_idx1][i] - spray_cloud2.uu[p_idx2][i];
      }
      rel_vel = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

      // decide which parcel has more drops
      CONVERGE_precision_t num_drop_max = spray_cloud1.num_drop[p_idx1];
      if(spray_cloud2.num_drop[p_idx2] > num_drop_max)
      {
	 num_drop_max = spray_cloud2.num_drop[p_idx2];
      }

      // calculate the chance of the two parcels being accepted for collision
      CONVERGE_precision_t prob_ij = (num_drop_max * rel_vel * cross_section) / max_vel_cross;
      random1                      = CONVERGE_random_precision();
      CONVERGE_precision_t random2 = CONVERGE_random_precision();

      if(random1 <= prob_ij)   // a collision occurs
      {
	 // set p_max to the parcel with the larger radius,
	 // set p_min to the parcel with the smaller radius
	 CONVERGE_precision_t p1_rad_save = spray_cloud1.radius[p_idx1];
	 CONVERGE_precision_t p2_rad_save = spray_cloud2.radius[p_idx2];

	 CONVERGE_int_t p_max, p_min, pc_max, pc_min;
	 if(p1_rad_save > p2_rad_save)
	 {
	    p_max  = p_idx1;
	    p_min  = p_idx2;
	    pc_max = pc_idx1;
	    pc_min = pc_idx2;
	 }
	 else
	 {
	    p_max  = p_idx2;
	    p_min  = p_idx1;
	    pc_max = pc_idx2;
	    pc_min = pc_idx1;
	 }

	 CONVERGE_cloud_t cvg_cloud_max = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_max);
	 struct ParcelCloud spray_cloud_max;
	 load_user_cloud_for_collision(&spray_cloud_max, cvg_cloud_max, parcel_type);

	 CONVERGE_cloud_t cvg_cloud_min = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_min);
	 struct ParcelCloud spray_cloud_min;
	 load_user_cloud_for_collision(&spray_cloud_min, cvg_cloud_min, parcel_type);

	 if(parcel_type == liquid_type)
	 {
	    CONVERGE_precision_t rad_ratio = spray_cloud_max.radius[p_max] / spray_cloud_min.radius[p_min];

	    CONVERGE_precision_t delt_post = 0.0, delt_post2 = 0.0, delt_post3 = 0.0, delt_post6 = 0.0;
	    if(flags->collision_outcome_flag == 1)
	    {
	       delt_post  = spray_cloud_min.radius[p_min] / spray_cloud_max.radius[p_max];
	       delt_post2 = delt_post * delt_post;
	       delt_post3 = delt_post2 * delt_post;
	       delt_post6 = delt_post3 * delt_post3;
	    }

	    CONVERGE_precision_t rad_func =
	       rad_ratio * rad_ratio * rad_ratio - 2.4 * rad_ratio * rad_ratio + 2.7 * rad_ratio;

	    // calculate mass averaged surface tension for the two drops
	    CONVERGE_precision_t rad_max3 = spray_cloud_max.density[p_max] * spray_cloud_max.radius[p_max] *
	       spray_cloud_max.radius[p_max] * spray_cloud_max.radius[p_max];
	    CONVERGE_precision_t rad_min3 = spray_cloud_min.density[p_min] * spray_cloud_min.radius[p_min] *
	       spray_cloud_min.radius[p_min] * spray_cloud_min.radius[p_min];

	    CONVERGE_precision_t surften_max = 0.0, surften_min = 0.0;
	    for(CONVERGE_int_t isp = 0; isp < num_parcel_species; isp++)
	    {
	       surften_max += spray_cloud_max.mfrac[p_max * num_parcel_species + isp] *
		  CONVERGE_table_lookup(surf_tension_table[isp], spray_cloud_max.temp[p_max]);

	       surften_min += spray_cloud_min.mfrac[p_min * num_parcel_species + isp] *
		  CONVERGE_table_lookup(surf_tension_table[isp], spray_cloud_min.temp[p_min]);
	    }

	    CONVERGE_precision_t surface_tension =
	       (rad_max3 * surften_max + rad_min3 * surften_min) / (rad_max3 + rad_min3);

	    // calculate collision weber number (note that Post's paper includes the sum of the two
	    // drop radii in the weber number, but we use the minimum radius as by O'Rourke)
	    CONVERGE_precision_t weber_collision =
	       spray_cloud_max.density[p_max] * rel_vel * rel_vel * spray_cloud_min.radius[p_min] / surface_tension;

	    // pterm is used to calculate b_crit below
	    CONVERGE_precision_t pterm = 2.4 * rad_func / (weber_collision + 1.0e-18);

	    // calculate value of b_crit^2/(r_1+r_2)^2 - Eq. 5.14 of O'Rourke's thesis
	    // this is E_coal in Post's paper
	    CONVERGE_precision_t prob_coal = (1.0 < pterm) ? 1.0 : pterm;

	    if(flags->collision_outcome_flag == 0)
	    {
	       random2 = CONVERGE_random_precision();
	       if(random2 <= prob_coal)
	       {
		  // coalescence occurs between the drops
		  CONVERGE_precision_t dummy1 = 1.0, dummy2 = 1.0, dummy3 = 1.0;

		  spray_cloud_min.from_injector =
		     (CONVERGE_int_t *)CONVERGE_cloud_get_field_data(cvg_cloud_min, PHASE_FROM_INJECTOR);
		  CONVERGE_injector_t injector = CONVERGE_get_injector_with_id(spray_cloud_min.from_injector[p_min]);
		  CONVERGE_precision_t inj_mass_per_parcel =
		     CONVERGE_injector_get_parameter_precision(injector, INJECTOR_MASS_PER_PARCEL);

		  CONVERGE_spray_coalesce(p_max,
					  p_min,
					  cvg_cloud_max,
					  cvg_cloud_min,
					  dummy1,
					  dummy2,
					  dummy3,
					  multflag,
					  inj_mass_per_parcel,
					  num_parcel_species);
	       }
	       else
	       {
		  // Grazing collision occurs
		  CONVERGE_precision_t dummy1 = 1.0;
		  CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, dummy1);
	       }
	    }
	    else
	    {
	       // calculate b_crit (impact_crit) value from prob_coal. also calculate b (impact_param)
	       // Post has the sum of the radii inside the sqrt for impact_param, but O'Rourke does not
	       CONVERGE_precision_t impact_crit =
		  sqrt(prob_coal) * (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);
	       CONVERGE_precision_t impact_param =
		  sqrt(random2) * (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);

	       // calculate bouncing criteria of Estrade et al. (brand_post is B in Post's paper)
	       CONVERGE_precision_t brand_post =
		  impact_param / (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);
	       CONVERGE_precision_t tao_post = (1.0 - brand_post) * (1.0 + delt_post);

	       CONVERGE_precision_t chi1 = 0.0;
	       if(tao_post < 1.0)
	       {
		  chi1 = 0.25 * tao_post * tao_post * (3.0 - tao_post);
	       }
	       else
	       {
		  chi1 = 1.0 - 0.25 * (2.0 - tao_post) * (2.0 - tao_post) * (1.0 + tao_post);
	       }

	       if(flags->collision_mesh_flag == 1)
	       {
		  cell_den = 0.5 * (spray_cloud1.gas_density[p_idx1] + spray_cloud2.gas_density[p_idx2]);
	       }

	       CONVERGE_precision_t phi_shape = 3.351 * CONVERGE_pow_2over3(cell_den / 1.16);

	       // 1.0-brand_post*brand_post in denominator is equivalent to (cos(arcsin(brand_post)))^2 in Post's paper
	       CONVERGE_precision_t webounce = (delt_post * (1.0 + delt_post2) * (4.0 * phi_shape - 12.0)) /
		  (chi1 * (1.0 - brand_post * brand_post) + 1.0e-10);

	       // determine which outcome will take place
	       // should this be weber based on radius?
	       CONVERGE_precision_t weberd = 2.0 * weber_collision;

	       if(weberd < webounce)
	       {
		  CONVERGE_precision_t fe_post = 0.0;
		  CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, fe_post);


		  continue;
	       }
	       else
	       {
		  // determine if stretching separation takes place
		  CONVERGE_precision_t bs_post =
		     impact_crit / (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);
		  CONVERGE_int_t istretch = 0;
		  if(brand_post > bs_post)
		  {
		     istretch = 1;
		  }

		  // reflexive separation of Ashgriz and Poo
		  CONVERGE_int_t ireflex       = 0;
		  CONVERGE_precision_t xi_post = 0.5 * brand_post * (1.0 + delt_post);
		  CONVERGE_precision_t eta1_post =
		     2.0 * (1.0 - xi_post) * (1.0 - xi_post) * sqrt(1.0 - xi_post * xi_post) - 1.0;
		  CONVERGE_precision_t eta2_post = -1.0;

		  // protect against sqrt of negative number
		  if(delt_post > xi_post)
		  {
		     eta2_post =
			2.0 * (delt_post - xi_post) * (delt_post - xi_post) * sqrt(delt_post2 - xi_post * xi_post) -
			delt_post3;
		  }

		  CONVERGE_precision_t web_ref;
		  if(eta1_post > 0.0 && eta2_post > 0.0)
		  {
		     CONVERGE_precision_t a_post = 7.0 * CONVERGE_pow_2over3(1.0 + delt_post3);
		     CONVERGE_precision_t b_post = 4.0 * (1.0 + delt_post2);
		     CONVERGE_precision_t c_post = delt_post * (1.0 + delt_post3) * (1.0 + delt_post3);
		     CONVERGE_precision_t d_post = delt_post6 * eta1_post + eta2_post;

		     web_ref = 3.0 * (a_post - b_post) * (c_post / d_post);

		     if(weberd > web_ref)
		     {
			ireflex = 1;
		     }
		  }

		  if(istretch == 1 && ireflex == 0)
		  {
		     CONVERGE_precision_t aa_post = (brand_post - bs_post) / (1.0 - bs_post);
		     CONVERGE_precision_t fe_post = 1.0 - aa_post * aa_post;

		     CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, fe_post);


		     continue;
		  }

		  if(istretch == 0 && ireflex == 1)
		  {
		     CONVERGE_precision_t fe_post = web_ref / (weberd + DBL_EPSILON);

		     CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, fe_post);


		     continue;
		  }

		  if(istretch == 1 && ireflex == 1)
		  {
		     CONVERGE_precision_t aa_post    = (brand_post - bs_post) / (1.0 - bs_post);
		     CONVERGE_precision_t fe_stretch = 1.0 - aa_post * aa_post;

		     CONVERGE_precision_t fe_ref = web_ref / (weberd + DBL_EPSILON);

		     CONVERGE_precision_t fe_post;
		     if(fe_stretch < fe_ref)
		     {
			fe_post = fe_stretch;
		     }
		     else
		     {
			fe_post = fe_ref;
		     }

		     CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, fe_post);


		     continue;
		  }

		  if(istretch == 0 && ireflex == 0)
		  {
		     // coalescence occurs between the drops
		     multflag                    = 0.0;
		     CONVERGE_precision_t dummy1 = 1.0, dummy2 = 1.0, dummy3 = 1.0;

		     spray_cloud_min.from_injector =
			(CONVERGE_int_t *)CONVERGE_cloud_get_field_data(cvg_cloud_min, PHASE_FROM_INJECTOR);
		     CONVERGE_injector_t injector = CONVERGE_get_injector_with_id(spray_cloud_min.from_injector[p_min]);
		     CONVERGE_precision_t inj_mass_per_parcel =
			CONVERGE_injector_get_parameter_precision(injector, INJECTOR_MASS_PER_PARCEL);

		     CONVERGE_spray_coalesce(p_max,
					     p_min,
					     cvg_cloud_max,
					     cvg_cloud_min,
					     dummy1,
					     dummy2,
					     dummy3,
					     multflag,
					     inj_mass_per_parcel,
					     num_parcel_species);
		  }
	       }
	    }
	 }
	 else
	 {
	    solid_parcel_outcome(cvg_cloud_min,
				 cvg_cloud_max,
				 &spray_cloud_min,
				 &spray_cloud_max,
				 p_min,
				 p_max,
				 random2,
				 num_parcel_species,
                                 flags);
	 }
      }
   }
}


static void parcel_collide_Orourke(CONVERGE_int_t cloud_size,
				   const CONVERGE_int_t *cloud_idx_array,
				   const CONVERGE_int_t *parcel_idx_array,
				   CONVERGE_precision_t cell_vol,
				   CONVERGE_precision_t cell_den,
				   CONVERGE_cloud_list_t spray_cloud_list,
				   CONVERGE_index_t num_parcel_species,
				   CONVERGE_precision_t dt,
				   CONVERGE_int_t parcel_type,
				   flags_t *flags)
{
   // loop over parcels in the current parcel cloud.
   // Compare each parcel with every other parcel in the current parcel cloud.
   // Start on 2nd parcel in list.
   for(CONVERGE_int_t ii = 1; ii < cloud_size; ii++)
   {
      const CONVERGE_int_t pc_idx1 = cloud_idx_array[ii];
      const CONVERGE_int_t p_idx1  = parcel_idx_array[ii];

      CONVERGE_cloud_t cvg_cloud1 = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_idx1);
      struct ParcelCloud spray_cloud1;
      load_user_cloud_for_collision(&spray_cloud1, cvg_cloud1, parcel_type);
      
      // skip parcel collison if the LISA breakup model is on and it is a parent parcel
      if(parcel_type == liquid_type)
      {
	 spray_cloud1.from_injector =
	    (CONVERGE_int_t *)CONVERGE_cloud_get_field_data(cvg_cloud1, PHASE_FROM_INJECTOR);
	 spray_cloud1.parent = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud1, LAGRANGIAN_PARENT);

	 CONVERGE_injector_t injector1 = CONVERGE_get_injector_with_id(spray_cloud1.from_injector[p_idx1]);
	 CONVERGE_index_t lisa_flag    = CONVERGE_injector_get_parameter_flag(injector1, INJECTOR_LISA_FLAG);
	 if(lisa_flag == 1 && spray_cloud1.parent[p_idx1] == 1)
	 {
	    // Destroy the cloud that is created above.
	    continue;
	 }
      }

      for(CONVERGE_int_t jj = 0; jj < ii; jj++)
      {
	 const CONVERGE_int_t pc_idx2 = cloud_idx_array[jj];
	 const CONVERGE_int_t p_idx2  = parcel_idx_array[jj];

	 CONVERGE_cloud_t cvg_cloud2 = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_idx2);
	 struct ParcelCloud spray_cloud2;
	 load_user_cloud_for_collision(&spray_cloud2, cvg_cloud2, parcel_type);

	 // skip parcel collison if the LISA breakup model is on and it is a parent parcel
	 if(parcel_type == liquid_type)
	 {
	    spray_cloud2.from_injector =
	       (CONVERGE_int_t *)CONVERGE_cloud_get_field_data(cvg_cloud2, PHASE_FROM_INJECTOR);
	    spray_cloud2.parent = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud2, LAGRANGIAN_PARENT);

	    CONVERGE_injector_t injector2 = CONVERGE_get_injector_with_id(spray_cloud2.from_injector[p_idx2]);
	    CONVERGE_int_t lisa_flag2     = CONVERGE_injector_get_parameter_flag(injector2, INJECTOR_LISA_FLAG);
	    if(lisa_flag2 == 1 && spray_cloud2.parent[p_idx2] == 1)
	    {
	       // Destroy the cloud that is created above.
	       continue;
	    }
	 }

	 CONVERGE_precision_t random1 = CONVERGE_random_precision();
	 CONVERGE_precision_t random2 = CONVERGE_random_precision();

	 // calculate the relative velocity between the two parcels
	 CONVERGE_precision_t parcel_rel_vel = sqrt((spray_cloud2.uu[p_idx2][0] - spray_cloud1.uu[p_idx1][0]) *
						    (spray_cloud2.uu[p_idx2][0] - spray_cloud1.uu[p_idx1][0]) +
						    (spray_cloud2.uu[p_idx2][1] - spray_cloud1.uu[p_idx1][1]) *
						    (spray_cloud2.uu[p_idx2][1] - spray_cloud1.uu[p_idx1][1]) +
						    (spray_cloud2.uu[p_idx2][2] - spray_cloud1.uu[p_idx1][2]) *
						    (spray_cloud2.uu[p_idx2][2] - spray_cloud1.uu[p_idx1][2]));

	 // calculate the sum of the two parcels' radii
	 CONVERGE_precision_t parcel_rad_sum = spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2];

	 // tot_coll is used in the calculation of tot_max_coll below
	 CONVERGE_precision_t tot_coll = spray_cloud1.num_drop[p_idx1] * spray_cloud2.num_drop[p_idx2] * PI *
	    parcel_rad_sum * parcel_rad_sum * parcel_rel_vel * dt / cell_vol;

	 // if tot_coll is less than or equal to 0.0, increment the second parcel
	 if(tot_coll <= 0.0)
	 {
	    // Destroy the cloud that is created above.
	    continue;
	 }

	 // set pmax to the parcel with the larger radius, set pmin to the parcel with the smaller radius
	 CONVERGE_precision_t p1_rad_save = spray_cloud1.radius[p_idx1];
	 CONVERGE_precision_t p2_rad_save = spray_cloud2.radius[p_idx2];

	 CONVERGE_int_t p_max, p_min, pc_max, pc_min;
	 if(p1_rad_save > p2_rad_save)
	 {
	    p_max  = p_idx1;
	    p_min  = p_idx2;
	    pc_max = pc_idx1;
	    pc_min = pc_idx2;
	 }
	 else
	 {
	    p_max  = p_idx2;
	    p_min  = p_idx1;
	    pc_max = pc_idx2;
	    pc_min = pc_idx1;
	 }

	 CONVERGE_cloud_t cvg_cloud_max = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_max);
	 struct ParcelCloud spray_cloud_max;
	 load_user_cloud_for_collision(&spray_cloud_max, cvg_cloud_max, parcel_type);

	 CONVERGE_cloud_t cvg_cloud_min = CONVERGE_cloud_list_get_cloud_at(spray_cloud_list, pc_min);
	 struct ParcelCloud spray_cloud_min;
	 load_user_cloud_for_collision(&spray_cloud_min, cvg_cloud_min, parcel_type);

	 // tot_max_coll is the mean value of the collison probability, pg. 135 of O'Rourke's thesis
	 CONVERGE_precision_t tot_max_coll = tot_coll / spray_cloud_max.num_drop[p_max];

	 // no_coll is the probability of no collisions
	 CONVERGE_precision_t no_coll = exp(-tot_max_coll);

	 if(random1 <= no_coll)
	 {
	    // No collisions take place.
	    // Destroy the clouds created and move tto next parcel.
	    continue;
	 }

	 if(parcel_type == liquid_type)
	 {
	    // collision takes place, calculate Brazier-Smith function rad_func (f in Post and Abraham's paper)
	    // rad_ratio is gamma and delt_post is Delta in Post and Abraham's paper
	    CONVERGE_precision_t rad_ratio = spray_cloud_max.radius[p_max] / spray_cloud_min.radius[p_min];

	    CONVERGE_precision_t delt_post = 0.0, delt_post2 = 0.0, delt_post3 = 0.0, delt_post6 = 0.0;
	    if(flags->collision_outcome_flag == 1)
	    {
	       delt_post  = spray_cloud_min.radius[p_min] / spray_cloud_max.radius[p_max];
	       delt_post2 = delt_post * delt_post;
	       delt_post3 = delt_post2 * delt_post;
	       delt_post6 = delt_post3 * delt_post3;
	    }

	    CONVERGE_precision_t rad_func =
	       rad_ratio * rad_ratio * rad_ratio - 2.4 * rad_ratio * rad_ratio + 2.7 * rad_ratio;

	    // calculate mass averaged surface tension for the two drops
	    CONVERGE_precision_t rad_max3 = spray_cloud_max.density[p_max] * spray_cloud_max.radius[p_max] *
	       spray_cloud_max.radius[p_max] * spray_cloud_max.radius[p_max];

	    CONVERGE_precision_t rad_min3 = spray_cloud_min.density[p_min] * spray_cloud_min.radius[p_min] *
	       spray_cloud_min.radius[p_min] * spray_cloud_min.radius[p_min];

	    CONVERGE_precision_t surften_max = 0.0, surften_min = 0.0;
	    for(CONVERGE_int_t isp = 0; isp < num_parcel_species; isp++)
	    {
	       surften_max += spray_cloud_max.mfrac[p_max * num_parcel_species + isp] *
		  CONVERGE_table_lookup(surf_tension_table[isp], spray_cloud_max.temp[p_max]);

	       surften_min += spray_cloud_min.mfrac[p_min * num_parcel_species + isp] *
		  CONVERGE_table_lookup(surf_tension_table[isp], spray_cloud_min.temp[p_min]);
	    }

	    CONVERGE_precision_t surface_tension =
	       (rad_max3 * surften_max + rad_min3 * surften_min) / (rad_max3 + rad_min3);

	    // calculate collision weber number (note that Post's paper includes the sum of the two
	    // drop radii in the weber number, but we use the minimum radius as by O'Rourke)
	    CONVERGE_precision_t weber_collision = spray_cloud_max.density[p_max] * parcel_rel_vel * parcel_rel_vel *
	       spray_cloud_min.radius[p_min] / surface_tension;

	    // pterm is used to calculate b_crit below
	    CONVERGE_precision_t pterm = 2.4 * rad_func / (weber_collision + 1.0e-18);

	    // calculate value of b_crit^2/(r_1+r_2)^2 - Eq. 5.14 of O'Rourke's thesis
	    // this is E_coal in Post's paper
	    CONVERGE_precision_t prob_coal = (1.0 < pterm) ? 1.0 : pterm;

	    if(flags->collision_outcome_flag == 0)
	    {
	       if(random2 <= prob_coal)
	       {
		  // coalescence occurs between the drops
		  CONVERGE_precision_t multflag = 1.0;

		  spray_cloud_min.from_injector =
		     (CONVERGE_int_t *)CONVERGE_cloud_get_field_data(cvg_cloud_min, PHASE_FROM_INJECTOR);
		  CONVERGE_injector_t injector_min =
		     CONVERGE_get_injector_with_id(spray_cloud_min.from_injector[p_min]);
		  CONVERGE_precision_t inj_mass_per_parcel =
		     CONVERGE_injector_get_parameter_precision(injector_min, INJECTOR_MASS_PER_PARCEL);

		  CONVERGE_spray_coalesce(p_max,
					  p_min,
					  cvg_cloud_max,
					  cvg_cloud_min,
					  no_coll,
					  tot_max_coll,
					  random1,
					  multflag,
					  inj_mass_per_parcel,
					  num_parcel_species);
	       }
	       else
	       {
		  // grazing collision occurs
		  CONVERGE_precision_t dummy1 = 1.0;
		  CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, dummy1);
	       }
	    }
	    else
	    {
	       // calculate b_crit (impact_crit) value from prob_coal. also calculate b (impact_param)
	       // Post has the sum of the radii inside the sqrt for impact_param, but O'Rourke does not
	       CONVERGE_precision_t impact_crit =
		  sqrt(prob_coal) * (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);
	       CONVERGE_precision_t impact_param =
		  sqrt(random2) * (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);

	       // calculate bouncing criteria of Estrade et al. (brand_post is B in Post's paper)
	       CONVERGE_precision_t brand_post =
		  impact_param / (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);
	       CONVERGE_precision_t tao_post = (1.0 - brand_post) * (1.0 + delt_post);

	       CONVERGE_precision_t chi1 = 0.0;
	       if(tao_post < 1.0)
	       {
		  chi1 = 0.25 * tao_post * tao_post * (3.0 - tao_post);
	       }
	       else
	       {
		  chi1 = 1.0 - 0.25 * (2.0 - tao_post) * (2.0 - tao_post) * (1.0 + tao_post);
	       }

	       if(flags->collision_mesh_flag == 1)
	       {
		  cell_den = 0.5 * (spray_cloud1.gas_density[p_idx1] + spray_cloud2.gas_density[p_idx2]);
	       }

	       CONVERGE_precision_t phi_shape = 3.351 * CONVERGE_pow_2over3(cell_den / 1.16);

	       // 1.0-brand_post*brand_post in denominator is equivalent to (cos(arcsin(brand_post)))^2 in Post's paper
	       CONVERGE_precision_t webounce = (delt_post * (1.0 + delt_post2) * (4.0 * phi_shape - 12.0)) /
		  (chi1 * (1.0 - brand_post * brand_post) + 1.0e-10);

	       // determine which outcome will take place
	       // should this be weber based on radius?
	       CONVERGE_precision_t weberd = 2.0 * weber_collision;

	       if(weberd < webounce)
	       {
		  CONVERGE_precision_t fe_post = 0.0;
		  CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, fe_post);

		  // Destroy the clouds created and before moving to next parcel.

		  continue;
	       }
	       else
	       {
		  // determine if stretching separation takes place
		  CONVERGE_precision_t bs_post =
		     impact_crit / (spray_cloud1.radius[p_idx1] + spray_cloud2.radius[p_idx2]);
		  CONVERGE_int_t istretch = 0;
		  if(brand_post > bs_post)
		  {
		     istretch = 1;
		  }

		  // reflexive separation of Ashgriz and Poo
		  CONVERGE_int_t ireflex       = 0;
		  CONVERGE_precision_t xi_post = 0.5 * brand_post * (1.0 + delt_post);
		  CONVERGE_precision_t eta1_post =
		     2.0 * (1.0 - xi_post) * (1.0 - xi_post) * sqrt(1.0 - xi_post * xi_post) - 1.0;
		  CONVERGE_precision_t eta2_post = -1.0;

		  // protect against sqrt of negative number
		  if(delt_post > xi_post)
		  {
		     eta2_post =
			2.0 * (delt_post - xi_post) * (delt_post - xi_post) * sqrt(delt_post2 - xi_post * xi_post) -
			delt_post3;
		  }

		  CONVERGE_precision_t web_ref;
		  if(eta1_post > 0.0 && eta2_post > 0.0)
		  {
		     CONVERGE_precision_t a_post = 7.0 * CONVERGE_pow_2over3(1.0 + delt_post3);
		     CONVERGE_precision_t b_post = 4.0 * (1.0 + delt_post2);
		     CONVERGE_precision_t c_post = delt_post * (1.0 + delt_post3) * (1.0 + delt_post3);
		     CONVERGE_precision_t d_post = delt_post6 * eta1_post + eta2_post;

		     web_ref = 3.0 * (a_post - b_post) * (c_post / d_post);

		     if(weberd > web_ref)
		     {
			ireflex = 1;
		     }
		  }

		  if(istretch == 1 && ireflex == 0)
		  {
		     CONVERGE_precision_t aa_post = (brand_post - bs_post) / (1.0 - bs_post);
		     CONVERGE_precision_t fe_post = 1.0 - aa_post * aa_post;

		     CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, fe_post);

		     // Destroy the clouds created and before moving to next parcel.

		     continue;
		  }

		  if(istretch == 0 && ireflex == 1)
		  {
		     CONVERGE_precision_t fe_post = web_ref / (weberd + DBL_EPSILON);
		     CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, fe_post);

		     // Destroy the clouds created and before moving to next parcel.

		     continue;
		  }

		  if(istretch == 1 && ireflex == 1)
		  {
		     CONVERGE_precision_t aa_post    = (brand_post - bs_post) / (1.0 - bs_post);
		     CONVERGE_precision_t fe_stretch = 1.0 - aa_post * aa_post;

		     CONVERGE_precision_t fe_ref = web_ref / (weberd + DBL_EPSILON);

		     CONVERGE_precision_t fe_post;
		     if(fe_stretch < fe_ref)
		     {
			fe_post = fe_stretch;
		     }
		     else
		     {
			fe_post = fe_ref;
		     }
		     CONVERGE_spray_graze_collide(p_idx1, p_idx2, cvg_cloud1, cvg_cloud2, prob_coal, random2, fe_post);

		     // Destroy the clouds created and before moving to next parcel.

		     continue;
		  }

		  if(istretch == 0 && ireflex == 0)
		  {
		     // coalescence occurs between the drops
		     CONVERGE_precision_t multflag = 1.0;

		     spray_cloud_min.from_injector =
			(CONVERGE_int_t *)CONVERGE_cloud_get_field_data(cvg_cloud_min, PHASE_FROM_INJECTOR);
		     CONVERGE_injector_t injector_min =
			CONVERGE_get_injector_with_id(spray_cloud_min.from_injector[p_min]);
		     CONVERGE_precision_t inj_mass_per_parcel =
			CONVERGE_injector_get_parameter_precision(injector_min, INJECTOR_MASS_PER_PARCEL);

		     CONVERGE_spray_coalesce(p_max,
					     p_min,
					     cvg_cloud_max,
					     cvg_cloud_min,
					     no_coll,
					     tot_max_coll,
					     random1,
					     multflag,
					     inj_mass_per_parcel,
					     num_parcel_species);
		  }
	       }
	    }
	 }
	 else
	 {
	    solid_parcel_outcome(cvg_cloud_min,
				 cvg_cloud_max,
				 &spray_cloud_min,
				 &spray_cloud_max,
				 p_min,
				 p_max,
				 random2,
				 num_parcel_species,
                                 flags);
	 }
	 // Destroy the clouds created and before moving to next parcel.
      }
   }
}

static void load_user_cloud_for_collision(struct ParcelCloud *spray_cloud_loc, CONVERGE_cloud_t cvg_cloud, CONVERGE_int_t parcel_type)
{
   if(parcel_type == liquid_type)
   {
      spray_cloud_loc->num_drop = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_NUM_DROP);
      spray_cloud_loc->radius   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_RADIUS);
      spray_cloud_loc->uu       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_UU);
      spray_cloud_loc->temp     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_TEMP);
      spray_cloud_loc->density  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_DENSITY);
      spray_cloud_loc->gas_density =
	 (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_GAS_DENSITY);
      spray_cloud_loc->mfrac = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, LAGRANGIAN_MFRAC);
   }
   else
   {
      spray_cloud_loc->num_drop = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, SOLID_PARCEL_NUM_DROP);
      spray_cloud_loc->radius   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, SOLID_PARCEL_RADIUS);
      spray_cloud_loc->uu       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(cvg_cloud, SOLID_PARCEL_UU);
      spray_cloud_loc->xx       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(cvg_cloud, SOLID_PARCEL_XX);
      spray_cloud_loc->temp     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, SOLID_PARCEL_TEMP);
      spray_cloud_loc->density  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, SOLID_PARCEL_DENSITY);
      spray_cloud_loc->mfrac = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(cvg_cloud, SOLID_PARCEL_MFRAC);
   }
}

static void initialize_tables(CONVERGE_species_t species)
{
   CONVERGE_iterator_t parcel_species_it;
   CONVERGE_species_parcel_iterator_create(species, &parcel_species_it);

   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_SURFACE_TENSION_ID, &surf_tension_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
}

static void destroy_tables(CONVERGE_species_t species)
{
   CONVERGE_iterator_t parcel_species_it;
   CONVERGE_species_parcel_iterator_create(species, &parcel_species_it);

   unload_species_tables(parcel_species_it, &surf_tension_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
}
