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
#include <globals.h>
#include <Vb.h>

// Global variables
// Global velocity variables
CONVERGE_precision_t user_child_velocity_x =0.0;
CONVERGE_precision_t user_child_velocity_y =0.0;
CONVERGE_precision_t user_child_velocity_z =0.0;


// Profiling accumulators
static double prof_calcs = 0.0;
static double prof_loop = 0.0;
static double prof_child_parcel = 0.0;
static double prof_property_copy = 0.0;
static double prof_velocity_calc = 0.0;
static double prof_insert_cloud = 0.0;
static int last_cycle = -1;
// Function to print profiling information


void Breakup(struct ParcelCloud *old_parcel_cloud, CONVERGE_index_t p_idx, CONVERGE_cloud_t cloud)
{
    // UDF-level check to prevent creating child parcels from a parent with non-physical properties
    if (isnan(old_parcel_cloud->temp[p_idx]) || isinf(old_parcel_cloud->temp[p_idx]) || old_parcel_cloud->temp[p_idx] < 100.0) {
        CONVERGE_logger_warn("Breakup.c: Parent parcel %d (cloud %d) has invalid temperature (%.2f) at ncyc %ld. Skipping breakup for this parcel.",
            p_idx, old_parcel_cloud->cloud_index[p_idx], CONVERGE_ncyc(), old_parcel_cloud->temp[p_idx]);
        return;
    }
    if (isnan(old_parcel_cloud->uu[p_idx][0]) || isinf(old_parcel_cloud->uu[p_idx][0]) ||
        isnan(old_parcel_cloud->uu[p_idx][1]) || isinf(old_parcel_cloud->uu[p_idx][1]) ||
        isnan(old_parcel_cloud->uu[p_idx][2]) || isinf(old_parcel_cloud->uu[p_idx][2])) {
        CONVERGE_logger_warn("Breakup.c: Parent parcel %d (cloud %d) has invalid velocity at ncyc %ld. Skipping breakup for this parcel.",
            p_idx, old_parcel_cloud->cloud_index[p_idx], CONVERGE_ncyc());
        return;
    }

    // printf("\n Breakup Triggered, r_bubble = %2e, v_bubble = %2e, radius = %2e, breakup_flag = %i\n", old_parcel_cloud->r_bubble[p_idx], old_parcel_cloud->v_bubble[p_idx], old_parcel_cloud->radius[p_idx], old_parcel_cloud->thermal_breakup_flag[p_idx]);

    //Timing vars
    CONVERGE_precision_t t0 = CONVERGE_mpi_wtime();
    // Section 1: Initialization and validation
    
    // Check if cloud and parcel cloud are valid
    if (!cloud || !old_parcel_cloud) {
        printf("\nBreakup.c: Invalid cloud or parcel cloud pointer\n");
        CONVERGE_mpi_abort();
    }
    CONVERGE_vec3_t user_child_velocity[12];

    // Get cloud size and verify parcel index
    CONVERGE_index_t cloud_size = CONVERGE_cloud_size(cloud);
    // printf("\nBreakup.c: Cloud size = %d, p_idx = %d\n", cloud_size, p_idx);
    if (p_idx >= cloud_size) {
        printf("\nBreakup.c: Invalid parcel index %ld (cloud size = %ld)\n", p_idx, cloud_size);
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
        printf("\nBreakup.c: Parent velocity at p_idx %ld is NULL\n", p_idx);
        CONVERGE_mpi_abort();
    }
    // printf("\nBreakup.c: Parent velocity = %e %e %e\n", 
        //    old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2]);

    // Verify child_uu exists for this parcel
    if (!old_parcel_cloud->child_uu[p_idx]) {
        printf("\nBreakup.c: child_uu at p_idx %ld is NULL\n", p_idx);
        CONVERGE_mpi_abort();
    }
    if(old_parcel_cloud->thermal_breakup_flag[p_idx]==4){
    printf("\n ERROR, breakup routine being triggered on parcel with tbf = 4, tbf = %i",old_parcel_cloud->thermal_breakup_flag[p_idx]);
    CONVERGE_mpi_abort();
}
 CONVERGE_precision_t   old_r = old_parcel_cloud->radius[p_idx];
 CONVERGE_precision_t  old_nd = old_parcel_cloud->num_drop[p_idx];
    // printf("\nbreakup count %i",breakup_counter);
    if (old_parcel_cloud->radius[p_idx] < old_parcel_cloud->r_bubble[p_idx])
    {
        old_parcel_cloud->r_bubble[p_idx] = 0.95 * old_parcel_cloud->radius[p_idx];
    }

    // old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    // printf("running thermal breakup routine p_idx = %i, breakup count = %i \n",p_idx);
    CONVERGE_index_t num_child_parcels =12;
    CONVERGE_index_t N = num_child_parcels;

    // End of initialization section
    // printf("\nbreakup count %i",breakup_counter);
    if (old_parcel_cloud->radius[p_idx] < old_parcel_cloud->r_bubble[p_idx])
    {
        old_parcel_cloud->r_bubble[p_idx] = 0.95 * old_parcel_cloud->radius[p_idx];
    }
    // Section 2: Velocity calculations
    
    // Create velocity vectors for all child parcels
    // Get parent parcel's velocity - v = vx i + vy j + vz k
    CONVERGE_precision_t parent_vx, parent_vy, parent_vz;
    CONVERGE_vec3_t parent_velocity,parent_velocity_unit,parent_normal ;
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], &parent_velocity);
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], &parent_velocity_unit); // Parent velocity vector

            // --- Safety check: clamp parent velocity if unphysical ---
        CONVERGE_precision_t parent_vel_mag = sqrt(parent_velocity[0]*parent_velocity[0] +
                                    parent_velocity[1]*parent_velocity[1] +
                                    parent_velocity[2]*parent_velocity[2]);

        // Define a threshold based on your expected flow speeds
        // (e.g. 1e3 m/s is likely way too high for sprays)
        const CONVERGE_precision_t parent_vel_limit = 1.0e3;

        if (parent_vel_mag > parent_vel_limit) {
            printf("Breakup.c WARNING: Parent velocity too large (%e). "
                "Clamping to zero.\n", parent_vel_mag);

            for (int i = 0; i < 3; i++) {
                if(fabs(parent_velocity[i])>parent_vel_limit){
                    parent_velocity[i] = parent_vel_limit;
                }
            }
        }




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
        // printf("\nLarge rad vel ---- parent vel magnitude = %e, child_rad_vel = %e", parent_vmag, rad_vel);
        // printf("\n p_idx = %i",p_idx);
        // printf("\n r_bubble = %e",old_parcel_cloud->r_bubble[p_idx]);
        // printf("\n radius = %e",old_parcel_cloud->radius[p_idx]);
        // printf("\n v_bubble = %e",old_parcel_cloud->v_bubble[p_idx]);   
        // printf("\n thermal_breakup_flag = %i",old_parcel_cloud->thermal_breakup_flag[p_idx]);
        // printf("\n tbt = %i",old_parcel_cloud->tbt[p_idx]);
        // printf("\n pbt = %e",old_parcel_cloud->pbt[p_idx]);
        // printf("\n parent velocity = %e %e %e",old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2]);
        // // rad_vel = 0.0;           //Tried this to fix probllem with YNH3 field but it didn't work
    
    }
    else if(fabs(rad_vel)<1.0e-9){
        printf("\n rad_vel = %e",rad_vel);
        printf("\n p_idx = %li",p_idx);
        printf("\n r_bubble = %e",old_parcel_cloud->r_bubble[p_idx]);
        printf("\n radius = %e",old_parcel_cloud->radius[p_idx]);
        printf("\n v_bubble = %e",old_parcel_cloud->v_bubble[p_idx]);   
        printf("\n thermal_breakup_flag = %i",old_parcel_cloud->thermal_breakup_flag[p_idx]);
        printf("\n tbt = %i",old_parcel_cloud->tbt[p_idx]);
        printf("\n pbt = %i",old_parcel_cloud->pbt[p_idx]);
        CONVERGE_mpi_abort();
    }


////Start of old child velocity _________________
    // printf("rad _vel =  %e, vmag = %e",rad_vel,parent_vmag);
    CONVERGE_precision_t aa = 10.0; // Scale factor for velocity 
 
//perpendicular vector calculation
CONVERGE_vec3_t arbitrary;
do {
    arbitrary[0] = 2.0 * CONVERGE_random_precision() - 1.0;
    arbitrary[1] = 2.0 * CONVERGE_random_precision() - 1.0;
    arbitrary[2] = 2.0 * CONVERGE_random_precision() - 1.0;
    CONVERGE_vec3_normalize(arbitrary);
} while (fabs(CONVERGE_vec3_dot(arbitrary, parent_velocity_unit)) > 0.7);


CONVERGE_vec3_normalize(arbitrary);
CONVERGE_vec3_cross(parent_velocity_unit, arbitrary, &parent_normal);
CONVERGE_vec3_normalize(parent_normal);
CONVERGE_precision_t normal_length = CONVERGE_vec3_length(parent_normal);



// Debug check - should be very close to 1.0
if (fabs(normal_length - 1.0) > 1.0e-1) {
    printf("ERROR: parent_normal not properly normalized! Length = %e\n", normal_length);
    printf("parent_normal = [%e, %e, %e]\n", 
           parent_normal[0], parent_normal[1], parent_normal[2]);
           printf("\n parent_velocity_unit = %e %e %e\n ", parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);
           printf("\n arbritrary = %e %e %e\n ", arbitrary[0], arbitrary[1], arbitrary[2]);
    // CONVERGE_mpi_abort();
}
    // printf("\nparent_normal = %e %e %e\n", parent_normal[0], parent_normal[1], parent_normal[2]);

    
   //-----------------------------Calculate child parcel velocities------------------------------------------------

    // First child parcel will have radial velocity along normal
    CONVERGE_vec3_dup(parent_normal,&user_child_velocity[0]); // Set first child parcel's velocity to be along the normal
    CONVERGE_vec3_normalize(user_child_velocity[0]);
    CONVERGE_vec3_scale(user_child_velocity[0], rad_vel * aa);
    


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

    CONVERGE_vec3_dup(user_child_velocity[0], &a); // Previous child parcel's velocity 

    CONVERGE_vec3_cross(parent_velocity_unit, user_child_velocity[0],&b); 
    CONVERGE_vec3_dup(parent_velocity_unit, &c); // Parent velocity unit vector
    
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
    for (int jj = 1; jj < N; jj++) // For all the other parcels 2:N
    {
    CONVERGE_vec3_dup(user_child_velocity[jj-1], &a); // Previous child parcel's velocity 

    CONVERGE_vec3_cross(parent_velocity_unit, user_child_velocity[jj-1],&b); 
    CONVERGE_vec3_dup(parent_velocity_unit, &c); // Parent velocity unit vector
    
    CONVERGE_vec3_scale(a, cos_theta) ; //Term 1 
    CONVERGE_vec3_scale(b, sin_theta); // Term 2 
    CONVERGE_vec3_scale(c,CONVERGE_vec3_dot(parent_velocity_unit, user_child_velocity[jj-1])* (1- cos_theta)); //Term 3

    CONVERGE_vec3_add(a,b,&d);
    CONVERGE_vec3_add(d,c, &user_child_velocity[jj]); // Final child velocity vector
    CONVERGE_vec3_normalize(user_child_velocity[jj]);
    if(rad_vel * aa > 100.0){
    //  printf("|Vc| = %f",rad_vel * aa);   
    }
    CONVERGE_vec3_scale(user_child_velocity[jj],rad_vel * aa);
    
    for (int k = 0; k < 3; k++)
    user_child_velocity[jj][k] *= (1.0 + 0.1 * (CONVERGE_random_precision() - 0.5));
   
    } // end of jj loop

  //End of old child initilisation 


  

    //----------------------------Calculate post breakup radius and number of drops for each child parcel----------------------------

    // Radius and Num Drop
    //Pre-Brekaup radius compare
   // printf("\nradius = %e  r_drop_0 = %e r_therm =%e dgre_cycles = %i",old_parcel_cloud->radius[p_idx],old_parcel_cloud->r_drop_0[p_idx],old_parcel_cloud->r_therm[p_idx],old_parcel_cloud->dgre_cycle_count[p_idx]);
    CONVERGE_precision_t rad_denom, rad_term1, rad_term2, rad_term3, rad_term4, parent_radius,parent_nd;
    parent_radius = old_parcel_cloud->radius[p_idx];
    parent_nd= old_parcel_cloud->num_drop[p_idx];
    CONVERGE_precision_t r_bubble_cube = CONVERGE_cube(old_parcel_cloud->r_bubble[p_idx]);
    
    CONVERGE_precision_t denom = CONVERGE_cube(old_parcel_cloud->radius[p_idx]) - r_bubble_cube;
    if (fabs(denom) < 1.0e-20) {
        printf("\nBreakup.c: Error: Denominator in rad_denom calculation is close to zero. Aborting.\n");
        CONVERGE_mpi_abort();
    }

    //this is a comment 

    rad_denom = 0.5 * (1.0 / denom);
    rad_term1 = (CONVERGE_square(old_parcel_cloud->radius[p_idx]) + CONVERGE_square(old_parcel_cloud->r_bubble[p_idx]));
    rad_term2 = 3.0 * CONVERGE_square(old_parcel_cloud->v_bubble[p_idx]) * (r_bubble_cube - (r_bubble_cube*old_parcel_cloud->r_bubble[p_idx]) * (1 /parent_radius));
    if (fabs(old_parcel_cloud->surf_ten[p_idx]) < 1.0e-12) {
        printf("\nBreakup.c: Error: Surface tension is close to zero (%e) for parcel %li. Aborting.\n", old_parcel_cloud->surf_ten[p_idx], p_idx);
        CONVERGE_mpi_abort();
    }
    rad_term3 = old_parcel_cloud->density[p_idx] / (3.0 * old_parcel_cloud->surf_ten[p_idx]);

    //printf("\nrad term 3 = %e, den = %e, surten = %e", rad_term3, old_parcel_cloud->density[p_idx], old_parcel_cloud->radius[p_idx]);
    // printf("\nTERM 2 V_BUBBLE = %e R_BUBBLE = %e R_DROP = %e",old_parcel_cloud->v_bubble[p_idx],old_parcel_cloud->r_bubble[p_idx],old_parcel_cloud->radius[p_idx]);
    rad_term4 = CONVERGE_square(rad_vel) / 2.0;
CONVERGE_precision_t radius_denominator = 2.0 * rad_denom * rad_term1 + rad_term3 * (rad_term2 * rad_denom - rad_term4);
if (fabs(radius_denominator) < 1.0e-20) {
    printf("\nBreakup.c: Error: Denominator for calculated_radius is close to zero (%e) for parcel %li. Aborting.\n", radius_denominator, p_idx);
    CONVERGE_mpi_abort();
}
CONVERGE_precision_t calculated_radius = 1.0 / radius_denominator;
// calculated_radius = calculated_radius * 0.1; // Testing decimating the radius further to see if intensifies evap issues
    // old_parcel_cloud->radius_tm1[p_idx] = old_parcel_cloud->radius[p_idx];
    if (calculated_radius < 0.0)
    {
        printf("\nBreakup.c: Error: calculated_radius is negative. Aborting.\n");
        CONVERGE_mpi_abort();
    }
    // printf("\ntbt=%i",old_parcel_cloud->tbt[p_idx]);
    if (old_parcel_cloud->r_bubble[p_idx] > parent_radius)
    {
        printf("\nbubble radius larger than droplet's original radius");
    }
    if (calculated_radius > parent_radius )
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



    prof_calcs += CONVERGE_mpi_wtime() - t0;
    
    //--------- Testing Child Parcel Introduction ----------------//
 
    // Calculate number of child parcels
    CONVERGE_precision_t old_mass, new_mass;
    CONVERGE_index_t nnn;
    CONVERGE_precision_t growth_rate, wave_length, radius_equil;
    CONVERGE_precision_t new_parcel_num_drop, new_parcel_mass, new_radius;
    CONVERGE_vec3_t new_parcel_uu;

    new_radius = calculated_radius ;
    // new_radius = old_parcel_cloud->radius[p_idx];
    //Calculate new number of droplets to conserve mass 
    old_mass = old_parcel_cloud->num_drop[p_idx] * 1.3333 * PI * CONVERGE_cube(old_parcel_cloud->radius[p_idx]);
    new_mass = old_mass / num_child_parcels;
    new_parcel_num_drop = new_mass / (1.3333 * PI * CONVERGE_cube(new_radius));
    // new_parcel_num_drop = old_parcel_cloud->num_drop[p_idx];

    //Try cooling the parcel to saturation temp -2 to prevent excessive evap 
    // Just doing this manually for 2 bar for now T -> 252 K 
    // old_parcel_cloud->temp[p_idx] = 252.0;


    t0 = CONVERGE_mpi_wtime();
    growth_rate = 0.0;
    wave_length = 0.0;
    CONVERGE_index_t initial_cloud_size = CONVERGE_cloud_size(cloud);
    // printf("\nInitial cloud size = %i",initial_cloud_size);
    // printf("\nParent parcel radius = %e, num_drop = %e", old_parcel_cloud->radius[p_idx], old_parcel_cloud->num_drop[p_idx]);
    if(initial_cloud_size >0)
    {
        CONVERGE_precision_t nd_before_break = old_parcel_cloud->num_drop[p_idx];

                    //Zero parent drop's radius - this triggers CONVERGE to remove the parent 

                    old_parcel_cloud->radius[p_idx] = 0.0; 
                    old_parcel_cloud->radius_tm1[p_idx] = 0.0;
                    old_parcel_cloud->num_drop[p_idx] = 0.0; 
                    old_parcel_cloud->num_drop_tm1[p_idx] = 0.0; 
                    // old_parcel_cloud->xx[p_idx][0] = 1.0; // This put's the parcel outside of the domain, so it will be removed 
            for(nnn = 0; nnn < num_child_parcels; nnn++)
            {
              
                CONVERGE_vec3_add(parent_velocity, user_child_velocity[nnn], &new_parcel_uu);

                if (CONVERGE_vec3_length(new_parcel_uu) > 1.0e3) {
                    printf("\nBreakup.c: Child velocity too large = %e %e %e. Capping to parent velocity.", new_parcel_uu[0], new_parcel_uu[1], new_parcel_uu[2]);
                    CONVERGE_vec3_dup(parent_velocity, &new_parcel_uu);
                }

                int rank;
	            CONVERGE_mpi_comm_rank(&rank);
      
                user_child_velocity_x = new_parcel_uu[0];
                user_child_velocity_y = new_parcel_uu[1];
                user_child_velocity_z = new_parcel_uu[2];
        

                    CONVERGE_precision_t t1 = CONVERGE_mpi_wtime();
               CONVERGE_spray_child_parcel(new_parcel_uu,
                                            growth_rate,
                                            wave_length,
                                            new_radius, 
                                            new_parcel_num_drop,
                                            p_idx,
                                            cloud);
                prof_child_parcel += CONVERGE_mpi_wtime() - t1;

            // reload after adding parcels
            load_user_cloud(old_parcel_cloud, cloud);
            }


            CONVERGE_index_t new_cloud_size = CONVERGE_cloud_size(cloud);
            // printf("\nNew cloud size = %i\n\n",new_cloud_size);
            // if(new_cloud_size <= initial_cloud_size)
            // {
            //     printf("\nError: New cloud size is not larger than initial cloud size after breakup");
            //     CONVERGE_mpi_abort();
            // }
        }
            // --------- End of Testing Child Parcel Introduction ----------------//

            prof_loop += CONVERGE_mpi_wtime() - t0;

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
    // End of child parcel section
    



    int ncyc = CONVERGE_ncyc();
    if (ncyc != last_cycle) {
    int rank;
    CONVERGE_mpi_comm_rank(&rank);
    // printf("Rank %d, Cycle %d Breakup profiling (s): calc=%%f, child_parcel=%%f, loop=%%f\n",
        //    rank, ncyc, prof_loop, prof_child_parcel, prof_calcs);

    double total = prof_loop + prof_calcs;
    // printf("Rank %d, Cycle %d Breakup profiling (%%): calc=%%f, child_parcel=%%f, loop=%%f\n",
        //    rank, ncyc,
        //    100.0*prof_calcs/total,
        //    100.0*prof_child_parcel/total,
        //    100.0*prof_loop/total);

    prof_loop = prof_child_parcel = prof_calcs = 0.0;
    last_cycle = ncyc;
}
}