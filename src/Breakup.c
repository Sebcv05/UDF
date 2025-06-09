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
#include <CrossProduct.h>
#include <spray_break.h>
#include <PsatNH3.h>
#include<Vb.h>
void Breakup(struct ParcelCloud *old_parcel_cloud, CONVERGE_index_t p_idx,CONVERGE_cloud_t cloud)
{

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
    CONVERGE_vec3_t child_velocity[20];
    struct childv
    {
        CONVERGE_precision_t vx[50];
        CONVERGE_precision_t vy[50];
        CONVERGE_precision_t vz[50];
    } c;

    // Create velocity vectors for all child parcels
    // Get parent parcel's velocity - v = vx i + vy j + vz k
    CONVERGE_precision_t parent_vx, parent_vy, parent_vz;
    CONVERGE_vec3_t parent_velocity,parent_velocity_unit;
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], parent_velocity);
    CONVERGE_vec3_dup(old_parcel_cloud->uu[p_idx], parent_velocity_unit); // Parent velocity vector
    CONVERGE_vec3_normalize(parent_velocity_unit);      //Parent Unit velocity vector 

    // printf("\n Breakup.c\n uu[p_idx] = %e %e %e\nparent_velocity = %e %e %e\nparent_velocity_unit = %e %e %e\n", old_parcel_cloud->uu[p_idx][0], old_parcel_cloud->uu[p_idx][1], old_parcel_cloud->uu[p_idx][2],parent_velocity[0], parent_velocity[1], parent_velocity[2], parent_velocity_unit[0], parent_velocity_unit[1], parent_velocity_unit[2]);


    // CONVERGE_precision_t parent_vmag = CONVERGE_sqrt(CONVERGE_square(parent_vx) + CONVERGE_square(parent_vy) + CONVERGE_square(parent_vz));
    CONVERGE_precision_t parent_vmag = CONVERGE_vec3_length(parent_velocity);
    // Calculate magnitude of child parcel velocity
    CONVERGE_precision_t rad_vel = 3.0 * old_parcel_cloud->v_bubble[p_idx] * CONVERGE_square(old_parcel_cloud->r_bubble[p_idx]) * (old_parcel_cloud->radius[p_idx] - old_parcel_cloud->r_bubble[p_idx]) / (CONVERGE_cube(old_parcel_cloud->radius[p_idx]) - CONVERGE_cube(old_parcel_cloud->r_bubble[p_idx]));
    if (rad_vel > parent_vmag)
    {
        printf("\nLarge rad vel ---- parent vel magnitude = %e, child_rad_vel = %e", parent_vmag, rad_vel);
    }
    // printf("rad _vel =  %e, vmag = %e",rad_vel,parent_vmag);
    CONVERGE_precision_t aa = 1; // Scale factor for velocity
   



    // Make unit vector perpendicular to direction ( a.b = 0)
    CONVERGE_precision_t parent_nmag = CONVERGE_sqrt(2 * CONVERGE_square(parent_velocity_unit[2]) + CONVERGE_square(parent_velocity_unit[0] + parent_velocity_unit[1])) / parent_velocity_unit[2];
    CONVERGE_vec3_t parent_normal;
    parent_normal[0] = 1 / parent_nmag; // Normalized x component
    parent_normal[1] = 1 / parent_nmag; // Normalized y component
    parent_normal[2] = -(parent_velocity_unit[0] + parent_velocity_unit[1]) / CONVERGE_sqrt(2 * CONVERGE_square(parent_velocity_unit[2]) + CONVERGE_square(parent_velocity_unit[0] + parent_velocity_unit[1])); // Normalized z component
    
    printf("\nparent_normal = %e %e %e\n", parent_normal[0], parent_normal[1], parent_normal[2]);

    
   //-----------------------------Calculate child parcel velocities----------------------------

    // First child parcel will have radial velocity along normal
    CONVERGE_vec3_dup(parent_normal,child_velocity[0]); // Set first child parcel's velocity to be along the normal
    


    // Rotate first parcel's velocity by a random angle
    CONVERGE_precision_t rand = CONVERGE_random_precision(); // Generates a random number between 0 and 1
    CONVERGE_precision_t psi = rand * 2 * PI;                // Random angle between 0 and 2 PI
    // printf("psi = %f deg\n",psi*360/(2*PI));
    //  Each of the remaining child parcels will be evenly distributed around the plane perpendicular to the normal
    CONVERGE_precision_t sin_psi = sin(psi);
    CONVERGE_precision_t cos_psi = cos(psi);
    //Rotation around parent's velocity vector by angle psi - Rodrigues' rotation formula
    // child velocity[i] = child_velocity[i-1]*cos(psi) + parent_velocity_normal_x_child_velocity[i-1] * sin(psi) + parent_velocity_x_parent_velocity_x_child_velocity[i-1] * (1 - cos(psi))
    CONVERGE_vec3_t a,b,c,d;

    CONVERGE_vec3_dup(child_velocity[0], a); // Previous child parcel's velocity 

    CONVERGE_vec3_cross(parent_velocity_unit, child_velocity[0],&b); 
    CONVERGE_vec3_dup(parent_velocity_unit, c); // Parent velocity unit vector
    
    CONVERGE_vec3_scale(a, cos_psi) ; //Term 1 
    CONVERGE_vec3_scale(b, sin_psi); // Term 2 
    CONVERGE_vec3_scale(c,CONVERGE_vec3_dot(parent_velocity_unit, child_velocity[0])* (1- cos_psi)) //Term 3

    CONVERGE_vec3_add(a,b,&d);
    CONVERGE_vec3_add(d,c, &child_velocity[0]); // Final child velocity vector

    CONVERGE_precision_t sin_theta, cos_theta, theta;
    theta = 2 * PI / N; // Angle between child parcels
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    //Developed to split parent parcel into N smaller parcels at breakup 
    for (int jj = 1; jj < N; jj++) // For all the other parcels 2:N
    {
    CONVERGE_vec3_dup(child_velocity[jj-1], a); // Previous child parcel's velocity 

    CONVERGE_vec3_cross(parent_velocity_unit, child_velocity[jj-1],&b); 
    CONVERGE_vec3_dup(parent_velocity_unit, c); // Parent velocity unit vector
    
    CONVERGE_vec3_scale(a, cos_theta) ; //Term 1 
    CONVERGE_vec3_scale(b, sin_theta); // Term 2 
    CONVERGE_vec3_scale(c,CONVERGE_vec3_dot(parent_velocity_unit, child_velocity[jj-1])* (1- cos_theta)) //Term 3

    CONVERGE_vec3_add(a,b,&d);
    CONVERGE_vec3_add(d,c, &child_velocity[jj]); // Final child velocity vector
   
   
   
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
    printf("\n Testing Child Parcel Introduction....\n");
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
                //Calcualte velocity of each child parcel

                CONVERGE_vec3_dup(child_velocity[nnn], new_parcel_uu); // Copy child's velocity from child_velocity array

                // CONVERGE_vec3_dup(new_parcel_uu, old_parcel_cloud->child_uu[p_idx][nnn]); // Store child's velocity in old_parcel_cloud


                // old_parcel_cloud->child_uu[p_idx][0] = c.vx[nnn]; // Store child's velocity direction so child can be displaced
                // old_parcel_cloud->child_uu[p_idx][1] = c.vy[nnn]; // Store child's velocity direction so child can be displaced
                // old_parcel_cloud->child_uu[p_idx][2] = c.vz[nnn]; // Store child's velocity direction so child can be displaced

               CONVERGE_spray_child_parcel(new_parcel_uu,
                                           growth_rate,
                                           wave_length,
                                           new_radius,
                                           new_parcel_num_drop,
                                           p_idx,
                                           cloud);
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