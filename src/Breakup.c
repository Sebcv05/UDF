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

    breakup_counter++;
    // printf("\nbreakup count %i",breakup_counter);
    if (old_parcel_cloud->radius[p_idx] < old_parcel_cloud->r_bubble[p_idx])
    {
        old_parcel_cloud->r_bubble[p_idx] = 0.95 * old_parcel_cloud->radius[p_idx];
    }

    // old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    // printf("running thermal breakup routine p_idx = %i, breakup count = %i \n",p_idx);
    CONVERGE_index_t N = 12;
    CONVERGE_precision_t N_pres = 12; // N in precision form
    // Calculate velocity
    struct childv
    {
        CONVERGE_precision_t vx[50];
        CONVERGE_precision_t vy[50];
        CONVERGE_precision_t vz[50];
    } c;

    // Create velocity vectors for all child parcels
    // Get parent parcel's velocity - v = vx i + vy j + vz k
    CONVERGE_precision_t parent_vx, parent_vy, parent_vz;
    parent_vx = old_parcel_cloud->uu[p_idx][0];
    parent_vy = old_parcel_cloud->uu[p_idx][1];
    parent_vz = old_parcel_cloud->uu[p_idx][2];

    CONVERGE_precision_t parent_vmag = CONVERGE_sqrt(CONVERGE_square(parent_vx) + CONVERGE_square(parent_vy) + CONVERGE_square(parent_vz));
    // Calculate magnitude of child parcel velocity
    CONVERGE_precision_t rad_vel = 3.0 * old_parcel_cloud->v_bubble[p_idx] * CONVERGE_square(old_parcel_cloud->r_bubble[p_idx]) * (old_parcel_cloud->radius[p_idx] - old_parcel_cloud->r_bubble[p_idx]) / (CONVERGE_cube(old_parcel_cloud->radius[p_idx]) - CONVERGE_cube(old_parcel_cloud->r_bubble[p_idx]));
    if (rad_vel > parent_vmag)
    {
        printf("\nLarge rad vel ---- parent vel magnitude = %e, child_rad_vel = %e", parent_vmag, rad_vel);
    }
    // printf("rad _vel =  %e, vmag = %e",rad_vel,parent_vmag);
    CONVERGE_precision_t aa = 1; // Scale factor for velocity
    // Make unit vector for direction
    CONVERGE_precision_t parent_vxu, parent_vyu, parent_vzu;
    parent_vxu = parent_vx / parent_vmag;
    parent_vyu = parent_vy / parent_vmag;
    parent_vzu = parent_vz / parent_vmag;

    // Make unit vector perpendicular to direction ( a.b = 0)
    CONVERGE_precision_t parent_nmag = CONVERGE_sqrt(2 * CONVERGE_square(parent_vzu) + CONVERGE_square(parent_vxu + parent_vyu)) / parent_vzu;
    CONVERGE_precision_t parent_nxu, parent_nyu, parent_nzu;
    parent_nxu = 1 / parent_nmag;
    parent_nyu = 1 / parent_nmag;
    parent_nzu = -(parent_vxu + parent_vyu) / CONVERGE_sqrt(2 * CONVERGE_square(parent_vzu) + CONVERGE_square(parent_vxu + parent_vyu));

    // First child parcel will have radial velocity along normal
    c.vx[0] = parent_nxu;
    c.vy[0] = parent_nyu;
    c.vz[0] = parent_nzu;

    struct vect k;
    struct vect v;
    struct vect k_x_v;
    struct vect k_x_k_x_v;
    // k is the axis of rotation - the direction of the parent parcel
    k.x = parent_vxu;
    k.y = parent_vyu;
    k.z = parent_vzu;
    // v is the vector to be rotated - the first parcel's new velocity
    v.x = parent_nxu;
    v.y = parent_nyu;
    v.z = parent_nzu;
    // Calculate k x v and k x (k x v)
    k_x_v = cross_product(k, v);
    k_x_k_x_v = cross_product(k, k_x_v);

    // Rotate first parcel's velocity by a random angle
    CONVERGE_precision_t rand = CONVERGE_random_precision(); // Generates a random number between 0 and 1
    CONVERGE_precision_t psi = rand * 2 * PI;                // Random angle between 0 and 2 PI
    // printf("psi = %f deg\n",psi*360/(2*PI));
    //  Each of the remaining child parcels will be evenly distributed around the plane perpendicular to the normal
    CONVERGE_precision_t sin_psi = sin(psi);
    CONVERGE_precision_t cos_psi = cos(psi);
    c.vx[0] = v.x + (k_x_v.x * sin_psi) + k_x_k_x_v.x * (1 - cos_psi);
    c.vy[0] = v.y + (k_x_v.y * sin_psi) + k_x_k_x_v.y * (1 - cos_psi);
    c.vz[0] = v.z + (k_x_v.z * sin_psi) + k_x_k_x_v.z * (1 - cos_psi);
    //Developed to split parent parcel into N smaller parcels at breakup 
    // for (int jj = 1; jj < N; jj++) // For all the other parcels 2:N
    // {
    //     // Using Rodrigues' rotation formula to rotate about parent velocity by angle theta
    //     // v_rot = v + (1-cos(theta)) k x (k x v) + sin(theta) k x v;
    //     CONVERGE_precision_t theta = jj * 2 * PI / N;
    //     // printf(" jj= %i, theta = %f\n",jj,(psi+theta) * 360 / (2*PI));
    //     //  Each of the remaining child parcels will be evenly distributed around the plane perpendicular to the normal
    //     c.vx[jj] = v.x + (k_x_v.x * sin(theta + psi)) + k_x_k_x_v.x * (1 - cos(theta + psi));
    //     c.vy[jj] = v.y + (k_x_v.y * sin(theta + psi)) + k_x_k_x_v.y * (1 - cos(theta + psi));
    //     c.vz[jj] = v.z + (k_x_v.z * sin(theta + psi)) + k_x_k_x_v.z * (1 - cos(theta + psi));
    // } // end of jj loop

    CONVERGE_vec3_t xx;
    CONVERGE_vec3_t uu;
    for (int kk = 0; kk < 3; kk++)
    {
        xx[kk] = old_parcel_cloud->xx[p_idx][kk];
    };
    // include a for loop to do it N times per parcel

    old_parcel_cloud->uu[p_idx][0] = c.vx[0] * aa * rad_vel + parent_vx;
    old_parcel_cloud->uu[p_idx][1] = c.vy[0] * aa * rad_vel + parent_vy;
    old_parcel_cloud->uu[p_idx][2] = c.vz[0] * aa * rad_vel + parent_vz;
    // uu = old_parcel_cloud->uu[p_idx];
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
    old_parcel_cloud->radius[p_idx] = 1.0 / (2.0 * rad_denom * rad_term1 + rad_term3 * (rad_term2 * rad_denom - rad_term4));
    old_parcel_cloud->radius_tm1[p_idx] = old_parcel_cloud->radius[p_idx];
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
    // Calculate number of child parcels
    CONVERGE_index_t num_child_parcels = 1;
    CONVERGE_index_t nnn;
    CONVERGE_precision_t growth_rate, wave_length, radius_equil;
    CONVERGE_precision_t new_parcel_num_drop, new_parcel_mass, new_radius;
    growth_rate = 0.0;
    wave_length = 0.0;
    CONVERGE_index_t initial_cloud_size = CONVERGE_cloud_size(cloud);
    printf("\nInitial cloud size = %i",initial_cloud_size);
    if(initial_cloud_size >0)
    {
            // for(nnn = 0; nnn < num_child_parcels; nnn++)
            // {
            //    CONVERGE_spray_child_parcel(old_parcel_cloud->uu[p_idx],
            //                                growth_rate,
            //                                wave_length,
            //                                0.1* old_parcel_cloud->radius[p_idx],
            //                                0.001 * old_parcel_cloud->num_drop[p_idx],
            //                                p_idx,
            //                                old_parcel_cloud);
            // }

            // reload after adding parcels
            // load_user_cloud(&old_parcel_cloud, old_parcel_cloud);
            // CONVERGE_index_t new_cloud_size = CONVERGE_cloud_size(old_parcel_cloud);
            // printf("\nNew cloud size = %i",new_cloud_size);
            if(new_cloud_size <= initial_cloud_size)
            {
                printf("\nError: New cloud size is not larger than initial cloud size after breakup");
                CONVERGE_mpi_abort();
            }
        }
            // --------- End of Testing Child Parcel Introduction ----------------//

        // if (old_parcel_cloud->radius[p_idx] < parent_radius )
        // {
        //     printf("radius decreased by thermal breakup ");
        // }

    // printf("\nN = %i",N);
    CONVERGE_precision_t old_nd;
    old_nd = old_parcel_cloud->num_drop[p_idx];
    old_parcel_cloud->num_drop[p_idx] = old_parcel_cloud->num_drop[p_idx] * CONVERGE_cube(parent_radius / old_parcel_cloud->radius[p_idx]);
    CONVERGE_precision_t mnew = old_parcel_cloud->num_drop[p_idx] * 1.3333 * PI * CONVERGE_cube(old_parcel_cloud->radius[p_idx]);
    //printf("\nm0 = %e m_old = %e m_new = %e",old_parcel_cloud->m0[p_idx],parent_nd*1.3333*PI*CONVERGE_cube(parent_radius),mnew);
    if(mnew > 1.01* old_parcel_cloud->m0[p_idx])
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