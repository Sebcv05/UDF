// Breakup.c

#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <Breakup.h>
#include <spray_break.h>
#include <PsatNH3.h>
#include<Vb.h>

static int user_velocity_index = 0;

void Breakup(struct ParcelCloud *old_parcel_cloud, CONVERGE_index_t p_idx,CONVERGE_cloud_t cloud)
{
    // Check if cloud and parcel cloud are valid
    if (!cloud || !old_parcel_cloud) {
        printf("\nBreakup.c: Invalid cloud or parcel cloud pointer\n");
        CONVERGE_mpi_abort();
    }
    CONVERGE_vec3_t user_child_velocity[20];

    // Get cloud size and verify parcel index
    CONVERGE_index_t cloud_size = CONVERGE_cloud_size(cloud);
    // printf("\nBreakup.c: Cloud size = %d, p_idx = %d\n", cloud_size, p_idx);
    if (p_idx >= cloud_size) {
        printf("\nBreakup.c: Invalid parcel index %d (cloud size = %d)\n", p_idx, cloud_size);
        CONVERGE_mpi_abort();
    }

    // Verify all required fields are loaded
    if (!old_parcel_cloud->child_uu) {
        printf("\nBreakup.c: child_uu field not loaded\n");
        CONVERGE_mpi_abort();
    }
    if (!old_parcel_cloud->uu) {
        printf("\nBreakup.c: uu field not loaded\n");
        CONVERGE_mpi_abort();
    }
    if (!old_parcel_cloud->radius) {
        printf("\nBreakup.c: radius field not loaded\n");
        CONVERGE_mpi_abort();
    }
    if (!old_parcel_cloud->r_bubble) {
        printf("\nBreakup.c: r_bubble field not loaded\n");
        CONVERGE_mpi_abort();
    }
    if (!old_parcel_cloud->v_bubble) {
        printf("\nBreakup.c: v_bubble field not loaded\n");
        CONVERGE_mpi_abort();
    }

    // Verify parent velocity exists
    if (!old_parcel_cloud->uu[p_idx]) {
        printf("\nBreakup.c: Parent velocity at p_idx %d is NULL\n", p_idx);
        CONVERGE_mpi_abort();
    }
    // printf("\nBreakup.c: Parent velocity = %e %e %e\n", 
        //    old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2]);

    // Verify child_uu exists for this parcel
    if (!old_parcel_cloud->child_uu[p_idx]) {
        printf("\nBreakup.c: child_uu at p_idx %d is NULL\n", p_idx);
        CONVERGE_mpi_abort();
    }
    if(old_parcel_cloud->thermal_breakup_flag[p_idx]==4){
    printf("\n ERROR, breakup routine being triggered on parcel with tbf = 4, tbf = %i",old_parcel_cloud->thermal_breakup_flag[p_idx]);
    CONVERGE_mpi_abort();
}
 CONVERGE_precision_t   old_r = old_parcel_cloud->radius[p_idx];
 CONVERGE_precision_t  old_nd = old_parcel_cloud->num_drop[p_idx];
    breakup_counter++;
    // printf("\nbreakup count %i",breakup_counter);
    if (old_parcel_cloud->radius[p_idx] < old_parcel_cloud->r_bubble[p_idx])
    {
        old_parcel_cloud->r_bubble[p_idx] = 0.95 * old_parcel_cloud->radius[p_idx];
    }

    // old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    // printf("running thermal breakup routine p_idx = %i, breakup count = %i \n",p_idx);
    CONVERGE_index_t N = 12;

    // Calculate velocity


    // Create velocity vectors for all child parcels
    // Get parent parcel's velocity - v = vx i + vy j + vz k
    CONVERGE_precision_t parent_vx, parent_vy, parent_vz;
    CONVERGE_vec3_t parent_velocity,parent_velocity_unit,parent_normal ;
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], parent_velocity);
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], parent_velocity_unit); // Parent velocity vector
    CONVERGE_vec3_normalize(parent_velocity_unit);      //Parent Unit velocity vector 
    if(CONVERGE_vec3_length(parent_velocity_unit)<0.99){
        printf("\n Breakup.c\n uu[p_idx] = %e %e %e\nparent_velocity = %e %e %e\nparent_velocity_unit = %e %e %e\n", old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2],parent_velocity[0], parent_velocity[1], parent_velocity[2], parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);
        printf("\n length of parent_velocity_unit = %e",CONVERGE_vec3_length(parent_velocity_unit));
        CONVERGE_mpi_abort();
    }else if(CONVERGE_vec3_length(parent_velocity_unit)>1.01){
        printf("\n Breakup.c\n uu[p_idx] = %e %e %e\nparent_velocity = %e %e %e\nparent_velocity_unit = %e %e %e\n", old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2],parent_velocity[0], parent_velocity[1], parent_velocity[2], parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);
        printf("\n length of parent_velocity_unit = %e",CONVERGE_vec3_length(parent_velocity_unit));
        CONVERGE_mpi_abort();
    }
    // printf("\n Breakup.c\n uu[p_idx] = %e %e %e\nparent_velocity = %e %e %e\nparent_velocity_unit = %e %e %e\n", old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2],parent_velocity[0], parent_velocity[1], parent_velocity[2], parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);


    // CONVERGE_precision_t parent_vmag = CONVERGE_sqrt(CONVERGE_square(parent_vx) + CONVERGE_square(parent_vy) + CONVERGE_square(parent_vz));
    CONVERGE_precision_t parent_vmag = CONVERGE_vec3_length(parent_velocity);
    // Calculate magnitude of child parcel velocity
    CONVERGE_precision_t rad_vel = 3.0 * old_parcel_cloud->v_bubble[p_idx] * CONVERGE_square(old_parcel_cloud->r_bubble[p_idx]) * (old_parcel_cloud->radius[p_idx] - old_parcel_cloud->r_bubble[p_idx]) / (CONVERGE_cube(old_parcel_cloud->radius[p_idx]) - CONVERGE_cube(old_parcel_cloud->r_bubble[p_idx]));
    if (rad_vel > parent_vmag)
    {
        printf("\nLarge rad vel ---- parent vel magnitude = %e, child_rad_vel = %e", parent_vmag, rad_vel);
        printf("\n p_idx = %i",p_idx);
        printf("\n r_bubble = %e",old_parcel_cloud->r_bubble[p_idx]);
        printf("\n radius = %e",old_parcel_cloud->radius[p_idx]);
        printf("\n v_bubble = %e",old_parcel_cloud->v_bubble[p_idx]);   
        printf("\n thermal_breakup_flag = %i",old_parcel_cloud->thermal_breakup_flag[p_idx]);
        printf("\n tbt = %i",old_parcel_cloud->tbt[p_idx]);
        printf("\n pbt = %e",old_parcel_cloud->pbt[p_idx]);
        printf("\n parent velocity = %e %e %e",old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2]);
        // rad_vel = 0.0;           //Tried this to fix probllem with YNH3 field but it didn't work
    
    }
    else if(fabs(rad_vel)<1.0e-9){
        printf("\n rad_vel = %e",rad_vel);
        printf("\n p_idx = %i",p_idx);
        printf("\n r_bubble = %e",old_parcel_cloud->r_bubble[p_idx]);
        printf("\n radius = %e",old_parcel_cloud->radius[p_idx]);
        printf("\n v_bubble = %e",old_parcel_cloud->v_bubble[p_idx]);   
        printf("\n thermal_breakup_flag = %i",old_parcel_cloud->thermal_breakup_flag[p_idx]);
        printf("\n tbt = %i",old_parcel_cloud->tbt[p_idx]);
        printf("\n pbt = %e",old_parcel_cloud->pbt[p_idx]);
        CONVERGE_mpi_abort();
    }
    // printf("rad _vel =  %e, vmag = %e",rad_vel,parent_vmag);
    CONVERGE_precision_t aa = 1; // Scale factor for velocity
   
// 
//perpendicular vector calculation
CONVERGE_vec3_t arbitrary;
if (fabs(parent_velocity_unit[0]) < 0.9) {
    arbitrary[0] = 1.0; arbitrary[1] = 0.0; arbitrary[2] = 0.0;
} else if (fabs(parent_velocity_unit[1]) < 0.9) {
    arbitrary[0] = 0.0; arbitrary[1] = 1.0; arbitrary[2] = 0.0;
} else {
    arbitrary[0] = 0.0; arbitrary[1] = 0.0; arbitrary[2] = 1.0;
}

// Pass address of parent_normal to store the cross product result
CONVERGE_vec3_cross(parent_velocity_unit, arbitrary, &parent_normal);

// Normalize the result
CONVERGE_precision_t normal_length = CONVERGE_vec3_normalize(parent_normal);

// Debug check - should be very close to 1.0
if (fabs(normal_length - 1.0) > 1.0e-2) {
    printf("ERROR: parent_normal not properly normalized! Length = %e\n", normal_length);
    printf("parent_normal = [%e, %e, %e]\n", 
           parent_normal[0], parent_normal[1], parent_normal[2]);
    CONVERGE_mpi_abort();
}
    // printf("\nparent_normal = %e %e %e\n", parent_normal[0], parent_normal[1], parent_normal[2]);

    
   //-----------------------------Calculate child parcel velocities------------------------------------------------

    // First child parcel will have radial velocity along normal
    CONVERGE_vec3_dup(parent_normal,&user_child_velocity[0]); // Set first child parcel's velocity to be along the normal
    


    // Rotate first parcel's velocity by a random angle
    CONVERGE_precision_t rand = CONVERGE_random_precision(); // Generates a random number between 0 and 1
    CONVERGE_precision_t psi = rand * 2 * PI;                // Random angle between 0 and 2 PI
    // printf("psi = %f deg\n",psi*360/(2*PI));
    //  Each of the remaining child parcels will be evenly distributed around the plane perpendicular to the normal
    CONVERGE_precision_t sin_psi = sin(psi);
    CONVERGE_precision_t cos_psi = cos(psi);
    //Rotation around parent's velocity vector by angle psi - Rodrigues' rotation formula
    // child velocity[i] = user_child_velocity[i-1]*cos(psi) + parent_velocity_normal_x_user_child_velocity[i-1] * sin(psi) + parent_velocity_x_parent_velocity_x_user_child_velocity[i-1] * (1 - cos(psi))
    CONVERGE_vec3_t a,b,c,d;

    CONVERGE_vec3_dup(user_child_velocity[0], a); // Previous child parcel's velocity 

    CONVERGE_vec3_cross(parent_velocity_unit, user_child_velocity[0],&b); 
    CONVERGE_vec3_dup(parent_velocity_unit, c); // Parent velocity unit vector
    
    CONVERGE_vec3_scale(a, cos_psi) ; //Term 1 
    CONVERGE_vec3_scale(b, sin_psi); // Term 2 
    CONVERGE_vec3_scale(c,CONVERGE_vec3_dot(parent_velocity_unit, user_child_velocity[0])* (1- cos_psi)); //Term 3

    CONVERGE_vec3_add(a,b,&d);
    CONVERGE_vec3_add(d,c, &user_child_velocity[0]); // Final child velocity vector

    CONVERGE_precision_t sin_theta, cos_theta, theta;
    theta = 2 * PI / N; // Angle between child parcels
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    //Developed to split parent parcel into N smaller parcels at breakup 
    for (int jj = 1; jj < N-1; jj++) // For all the other parcels 2:N
    {
    CONVERGE_vec3_dup(user_child_velocity[jj-1], a); // Previous child parcel's velocity 

    CONVERGE_vec3_cross(parent_velocity_unit, user_child_velocity[jj-1],&b); 
    CONVERGE_vec3_dup(parent_velocity_unit, c); // Parent velocity unit vector
    
    CONVERGE_vec3_scale(a, cos_theta) ; //Term 1 
    CONVERGE_vec3_scale(b, sin_theta); // Term 2 
    CONVERGE_vec3_scale(c,CONVERGE_vec3_dot(parent_velocity_unit, user_child_velocity[jj-1])* (1- cos_theta)); //Term 3

    CONVERGE_vec3_add(a,b,&d);
    CONVERGE_vec3_add(d,c, &user_child_velocity[jj]); // Final child velocity vector
   CONVERGE_vec3_normalize(user_child_velocity[jj]);
    CONVERGE_vec3_scale(user_child_velocity[jj],rad_vel * aa);
    
   
    } // end of jj loop



    //----------------------------Calculate post breakup radius and number of drops for each child parcel----------------------------

    // Radius and Num Drop
    //Pre-Brekaup radius compare
   // printf("\nradius = %e  r_drop_0 = %e r_therm =%e dgre_cycles = %i",old_parcel_cloud->radius[p_idx],old_parcel_cloud->r_drop_0[p_idx],old_parcel_cloud->r_therm[p_idx],old_parcel_cloud->dgre_cycle_count[p_idx]);
    CONVERGE_precision_t rad_denom, rad_term1, rad_term2, rad_term3, rad_term4, parent_radius,parent_nd;
    parent_radius = old_parcel_cloud->radius[p_idx];
    parent_nd= old_parcel_cloud->num_drop[p_idx];
    rad_denom = 0.5 * (1.0 / (CONVERGE_cube(old_parcel_cloud->radius[p_idx]) - CONVERGE_cube(old_parcel_cloud->r_bubble[p_idx])));
    rad_term1 = (CONVERGE_square(old_parcel_cloud->radius[p_idx]) + CONVERGE_square(old_parcel_cloud->r_bubble[p_idx]));
    rad_term2 = 3.0 * CONVERGE_square(old_parcel_cloud->v_bubble[p_idx]) * (CONVERGE_cube(old_parcel_cloud->r_bubble[p_idx]) - (old_parcel_cloud->r_bubble[p_idx] * old_parcel_cloud->r_bubble[p_idx] * old_parcel_cloud->r_bubble[p_idx] * old_parcel_cloud->r_bubble[p_idx]) * (1 / old_parcel_cloud->radius[p_idx]));
    rad_term3 = old_parcel_cloud->density[p_idx] / (3.0 * old_parcel_cloud->surf_ten[p_idx]);

    //printf("\nrad term 3 = %e, den = %e, surten = %e", rad_term3, old_parcel_cloud->density[p_idx], old_parcel_cloud->radius[p_idx]);
    // printf("\nTERM 2 V_BUBBLE = %e R_BUBBLE = %e R_DROP = %e",old_parcel_cloud->v_bubble[p_idx],old_parcel_cloud->r_bubble[p_idx],old_parcel_cloud->radius[p_idx]);
    rad_term4 = CONVERGE_square(rad_vel) / 2.0;
    // old_parcel_cloud->radius[p_idx] = 1.0 / (2.0 * rad_denom * rad_term1 + rad_term3 * (rad_term2 * rad_denom - rad_term4));
    // old_parcel_cloud->radius_tm1[p_idx] = old_parcel_cloud->radius[p_idx];
    if (old_parcel_cloud->radius[p_idx] < 0.0)
    {
        printf("\n radius negative \n");
    }
    // printf("\ntbt=%i",old_parcel_cloud->tbt[p_idx]);
    if (old_parcel_cloud->r_bubble[p_idx] > parent_radius)
    {
        printf("\nbubble radius larger than droplet's original radius");
    }
    if (old_parcel_cloud->radius[p_idx] > parent_radius )
    {
        CONVERGE_int_t rankb;
        CONVERGE_mpi_comm_rank(&rankb);
        printf("\n Thermal Breakup has increased radius");
        printf("\n rank = %i p_idx = %i",rankb,p_idx);
        printf("\nr_old = %e r_new = %e rb = %e vb = %e", parent_radius, old_parcel_cloud->radius[p_idx], old_parcel_cloud->r_bubble[p_idx], old_parcel_cloud->v_bubble[p_idx]);
       
        printf("\nrad_term1 = %e rt2 = %e rt3 = %e rt4 = %e rd = %e", rad_term1, rad_term2, rad_term3, rad_term4, rad_denom);
        printf("\ntbf = %i tbt = %i", old_parcel_cloud->thermal_breakup_flag[p_idx], old_parcel_cloud->tbt[p_idx]);
        printf("\n vb = %e ",old_parcel_cloud->v_bubble[p_idx]);
        printf("\n Vc = %e rho/sigma = %e", rad_vel, 3.0 * rad_term3);
        printf("\nrt1*rd = %e", rad_term1 * rad_denom);
        printf("\nrt2*rt3*rd = %e", rad_term2 * rad_term3 * rad_denom);
        printf("\nrt4*rt3 = %e", rad_term4 * rad_term3);
        printf("\n DGRE cycle count = %i",old_parcel_cloud->dgre_cycle_count[p_idx]);
        printf("\ntbt = %i", old_parcel_cloud->tbt[p_idx]);
        CONVERGE_precision_t P_sat;
        CONVERGE_precision_t Td = old_parcel_cloud->temp[p_idx];
        Saturation_PressureNH3(Td,&P_sat);
        printf("\n P_sat = %e,Td = %f",P_sat,Td);
         printf("\n recalculating v_bubble...");
         CONVERGE_precision_t P_amb = 1.5e6; 
         Bubble_Velocity(old_parcel_cloud,p_idx,P_sat,P_amb);
        printf("\n v_bubble = %e",old_parcel_cloud->v_bubble[p_idx]);
        // printf("\nsurf ten = %e, density = %e rt3 = %e",old_parcel_cloud->surf_ten[p_idx],old_parcel_cloud->density[p_idx],rad_term3);
    }

      
    //--------- Testing Child Parcel Introduction ----------------//
    // printf("\n Testing Child Parcel Introduction....\n");
//Check to see if reducing parent's temperature prevents convergence failure
    // old_parcel_cloud->temp[p_idx] = old_parcel_cloud->temp[p_idx] -10.0;

    // Calculate number of child parcels
    CONVERGE_precision_t old_mass, new_mass;
    CONVERGE_index_t nnn;
    CONVERGE_precision_t growth_rate, wave_length, radius_equil;
    CONVERGE_precision_t new_parcel_num_drop, new_parcel_mass, new_radius;
    CONVERGE_vec3_t new_parcel_uu;
    CONVERGE_index_t num_child_parcels = 10;
    // Child radius and number of drops
    new_radius = old_parcel_cloud->radius[p_idx] * 0.1;
    old_mass = old_parcel_cloud->num_drop[p_idx] * 1.3333 * PI * CONVERGE_cube(old_parcel_cloud->radius[p_idx]);
    new_mass = old_mass / num_child_parcels;
    new_parcel_num_drop = new_mass / (1.3333 * PI * CONVERGE_cube(new_radius));



    growth_rate = 0.0;
    wave_length = 0.0;
    CONVERGE_index_t initial_cloud_size = CONVERGE_cloud_size(cloud);
    // printf("\nInitial cloud size = %i",initial_cloud_size);
    // printf("\nParent parcel radius = %e, num_drop = %e", old_parcel_cloud->radius[p_idx], old_parcel_cloud->num_drop[p_idx]);
    if(initial_cloud_size >0)
    {
        CONVERGE_precision_t nd_before_break = old_parcel_cloud->num_drop[p_idx];
            for(nnn = 0; nnn < num_child_parcels; nnn++)
            {
                // Debug: Print memory addresses and values before operations
                // printf("\nBreakup.c: Before operations - p_idx = %d\n", p_idx);
                // printf("\nBreakup.c: parent_uu = %e %e %e\n ", old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2]);
                // printf("\nBreakup.c: parent_normal = %e %e %e\n", parent_normal[0], parent_normal[1], parent_normal[2]);
                // printf("\nBreakup.c: rad_vel = %e\n", rad_vel);
                // printf("\nBreakup.c: child_uu address = %p\n", (void*)&old_parcel_cloud->child_uu[p_idx]);
                // printf("\nBreakup.c: radius address = %p\n", (void*)&old_parcel_cloud->radius);
                // Store the radial velocity component in child_uu
                // The radial velocity is calculated as rad_vel * parent_normal
                // printf("\nBreakup.c: Storing radial velocity in child_uu\n");
                // CONVERGE_vec3_scale(parent_normal, rad_vel);  // This modifies parent_normal in place
                
                // Copy the modified parent_normal to child_uu
                // CONVERGE_vec3_dup(user_child_velocity[nnn],&old_parcel_cloud->child_uu[p_idx]);
                //Manually copy to child_uu since dup isn't working
             
                // printf("\nBreakup.c: user_child_velocity = %e %e %e\n", user_child_velocity[nnn][0], user_child_velocity[nnn][1], user_child_velocity[nnn][2]);
                // printf("\nBreakup.c: child_uu = %e %e %e\n", old_parcel_cloud->child_uu[p_idx][0], old_parcel_cloud->child_uu[p_idx][1], old_parcel_cloud->child_uu[p_idx][2]);
                // old_parcel_cloud->child_index[p_idx] = user_velocity_index;
                user_velocity_index = (user_velocity_index + 1) % 20;  // Wrap around if needed

                // Debug: Verify values after storing
                // printf("\nBreakup.c: After storing - child_uu = %e %e %e\n", 
                //        old_parcel_cloud->child_uu[p_idx][0], old_parcel_cloud->child_uu[p_idx][1], old_parcel_cloud->child_uu[p_idx][2]);
                
                // Calculate child's final velocity by adding parent velocity to radial component
                // printf("\nBreakup.c: Calculating final velocity\n");
                // CONVERGE_vec3_add(old_parcel_cloud->child_uu[p_idx], old_parcel_cloud->uu[p_idx], &new_parcel_uu);
                CONVERGE_vec3_add(old_parcel_cloud->uu[p_idx], user_child_velocity[nnn], &new_parcel_uu)
                // Debug: Verify final velocity
                // printf("\nBreakup.c: Final velocity = %e %e %e\n", 
                    //    new_parcel_uu[0], new_parcel_uu[1], new_parcel_uu[2]);
                
                // Debug print of radial velocity component stored in child_uu
                // printf("\nBreakup.c radial velocity = %e %e %e\n", old_parcel_cloud->child_uu[p_idx][0], old_parcel_cloud->child_uu[p_idx][1], old_parcel_cloud->child_uu[p_idx][2]);
                // printf("\nBreakup.c child velocity = %e %e %e\n", user_child_velocity[nnn][0], user_child_velocity[nnn][1], user_child_velocity[nnn][2]);
                // old_parcel_cloud->child_uu[p_idx][0] = c.vx[nnn]; // Store child's velocity direction so child can be displaced
                // old_parcel_cloud->child_uu[p_idx][1] = c.vy[nnn]; // Store child's velocity direction so child can be displaced
                // // old_parcel_cloud->child_uu[p_idx][2] = c.vz[nnn]; // Store child's velocity direction so child can be displaced
                // printf("\nBreakup.c: Before CONVERGE_spray_child_parcel - child_uu = %e %e %e at %p\n",
                //     old_parcel_cloud->child_uu[p_idx][0],
                //     old_parcel_cloud->child_uu[p_idx][1],
                //     old_parcel_cloud->child_uu[p_idx][2],
                //     (void*)&old_parcel_cloud->child_uu[p_idx]);


               CONVERGE_spray_child_parcel(new_parcel_uu,
                                            growth_rate,
                                            wave_length,
                                            new_radius,
                                            new_parcel_num_drop,
                                            p_idx,
                                            cloud);
                //Update position of child Parcel

                // CONVERGE_index_t new_cloud_size = CONVERGE_cloud_size(cloud);
                // printf("\nBreakup.c: initial_cloud_size = %i, new_cloud_size = %i\n", initial_cloud_size, new_cloud_size);
                // printf("\nBreakup.c: nnn = %i\n", nnn);
                // CONVERGE_index_t child_idx = initial_cloud_size+nnn;

                // if(child_idx < new_cloud_size) {
                //     CONVERGE_vec3_t child_displacement,old_position;
                //     CONVERGE_vec3_dup(user_child_velocity[nnn], &child_displacement);
                //     CONVERGE_vec3_normalize(child_displacement);
                //     CONVERGE_vec3_scale(child_displacement, parent_radius);
                //     CONVERGE_vec3_dup(old_parcel_cloud->xx[child_idx], &old_position);
                // CONVERGE_vec3_add(old_position, child_displacement, &old_parcel_cloud->xx[child_idx]);
                // }
            }

            // reload after adding parcels
            load_user_cloud(old_parcel_cloud, cloud);

            //Update parent drop's radius
            // old_parcel_cloud->radius[p_idx] = 0old_parcel_cloud->radius[p_idx];
            old_parcel_cloud->radius[p_idx] = new_radius; // Set parent radius to new radius
            old_parcel_cloud->num_drop[p_idx] = new_parcel_num_drop; // Set parent num_drop to new num_drop
            // old_parcel_cloud->radius[p_idx]=0.0;
            old_parcel_cloud->temp[p_idx] = 250.0;
            old_parcel_cloud->radius_tm1[p_idx] = old_parcel_cloud->radius[p_idx];
            old_parcel_cloud->num_drop[p_idx] = 0.1 * nd_before_break;
            old_parcel_cloud->num_drop[p_idx] = 0;
            old_parcel_cloud->pbt[p_idx] = 0;
            old_parcel_cloud->thermal_breakup_flag[p_idx] = 5; // Set to 5 to prevent secondary breakup
            old_parcel_cloud->tbt[p_idx] = 0; // Reset thermal breakup time
        
            CONVERGE_index_t new_cloud_size = CONVERGE_cloud_size(cloud);
            // printf("\nNew cloud size = %i\n\n",new_cloud_size);
            // if(new_cloud_size <= initial_cloud_size)
            // {
            //     printf("\nError: New cloud size is not larger than initial cloud size after breakup");
            //     CONVERGE_mpi_abort();
            // }
        }
            // --------- End of Testing Child Parcel Introduction ----------------//

        // if (old_parcel_cloud->radius[p_idx] < parent_radius )
        // {
        //     printf("radius decreased by thermal breakup ");
        // }

    // printf("\nN = %i",N);

    // old_parcel_cloud->num_drop[p_idx] = old_parcel_cloud->num_drop[p_idx] * CONVERGE_cube(parent_radius / old_parcel_cloud->radius[p_idx]);
    CONVERGE_precision_t mnew = old_parcel_cloud->num_drop[p_idx] * 1.3333 * PI * CONVERGE_cube(old_parcel_cloud->radius[p_idx]);
    //printf("\nm0 = %e m_old = %e m_new = %e",old_parcel_cloud->m0[p_idx],parent_nd*1.3333*PI*CONVERGE_cube(parent_radius),mnew);
    if(mnew > 0.1*1.01* old_parcel_cloud->m0[p_idx])
    {
        printf("\nBreakup Model has increased droplet mass!!!\n m_old = %e m_new = %e\nnd_old = %e nd_new = %e\nr_old = %e r_new = %e\n Aborting!!!!!!!!\n",old_parcel_cloud->m0[p_idx],mnew,old_nd,old_parcel_cloud->num_drop[p_idx],old_r,old_parcel_cloud->radius[p_idx]);
        CONVERGE_mpi_abort();
    }
    old_parcel_cloud->thermal_breakup_flag[p_idx] = 4;
    old_parcel_cloud->tbt[p_idx] = 0;
    // old_parcel_cloud->kb[p_idx]=0;
    old_parcel_cloud->int_omega[p_idx]=0.0;
    old_parcel_cloud->r_bubble[p_idx]=0.0;

    if (old_parcel_cloud->num_drop[p_idx] < 0.0)
    {
        printf("\nParcel N_drop < 0!!!!\n");
    }
    old_parcel_cloud->tbreak_kh[p_idx] = old_parcel_cloud->thermal_breakup_flag[p_idx];
}