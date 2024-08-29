// TABDistort.c
// This function calculates the TAB parameters y and ydot
#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <TABDistort.h>

/// @brief Calculates y and ydot based on an a modified density and viscosity for the bubble-droplet system
/// @param old_parcel_cloud - structure storing Lagrangian data
/// @param p_idx - parcel index in cloud
/// @param dt - timestep 
/// @param global_density - Gas Phase density at parcel's location
/// @param mu_v - Gas phase viscosity
/// @param rho_b - bubble density
/// @return - void
void TABDistort(struct ParcelCloud *old_parcel_cloud, CONVERGE_index_t p_idx, CONVERGE_precision_t dt, CONVERGE_precision_t global_density, CONVERGE_precision_t mu_v,CONVERGE_precision_t rho_b)
{
    
    // TAB  variables
    CONVERGE_precision_t spray_tab_csubd = CONVERGE_get_double("lagrangian.tab_csubd");
    CONVERGE_precision_t spray_tab_csubk = CONVERGE_get_double("lagrangian.tab_csubk");
    CONVERGE_precision_t spray_tab_cfocbck = CONVERGE_get_double("lagrangian.tab_cfocbck");
    CONVERGE_injector_t injector = CONVERGE_get_injector_with_id(old_parcel_cloud->from_injector[p_idx]);
    CONVERGE_precision_t rho_l = old_parcel_cloud->density[p_idx];
    CONVERGE_precision_t mu = old_parcel_cloud->viscosity[p_idx];

    /*if (CONVERGE_injector_get_tab_flag(injector))
    {
       continue;
    }

    if (CONVERGE_injector_get_lisa_flag(injector))
    {
       continue;
    }
  */
    CONVERGE_precision_t old_y = old_parcel_cloud->distort[p_idx];
    CONVERGE_precision_t old_ydot = old_parcel_cloud->distort_dot[p_idx];

    // Calculate void fraction if bubble exists
    CONVERGE_precision_t sys_density, sys_viscosity, void_fraction;
    if (old_parcel_cloud->r_bubble[p_idx] > 1.0e-20)
    {
        // printf("r_bubble>0\n");
        CONVERGE_precision_t inv_Delta = old_parcel_cloud->r_bubble[p_idx] / old_parcel_cloud->radius[p_idx];
        void_fraction = inv_Delta * inv_Delta * inv_Delta;
        sys_density = void_fraction * rho_b + (1 - void_fraction) * rho_l;
        sys_viscosity = void_fraction * mu_v + (1 - void_fraction) * mu;
        // printf("gas mu = %e  liquid mu = %e\n",mu_v,mu);
    }
    else if (old_parcel_cloud->r_bubble[p_idx] < 1.0e-20)
    {
        // printf("r_bubble = 0\n");
        void_fraction = 0.0;
        sys_density = rho_l;
        sys_viscosity = mu;
    }

    // calculate weber number
    CONVERGE_precision_t weber = global_density * old_parcel_cloud->rel_vel_mag[p_idx] *
                                 old_parcel_cloud->rel_vel_mag[p_idx] * old_parcel_cloud->radius[p_idx] /
                                 old_parcel_cloud->surf_ten[p_idx];
    // We_c
    CONVERGE_precision_t We_c = weber * spray_tab_cfocbck;
    // printf("\n We_c = %f, v_rel = %e, den = %e, rad = %e, st = %e",weber * spray_tab_cfocbck,old_parcel_cloud.rel_vel_mag[p_idx],global_density[node_index],old_parcel_cloud.radius[p_idx],old_parcel_cloud.surf_ten[p_idx]);

    // tab_rtd is 1/t_d in TAB paper referenced in header -- sys_viscosity used in place of mu_l, sys_density used in place of rho_l
    CONVERGE_precision_t tab_rtd =
        0.5 * spray_tab_csubd * sys_viscosity /
        (sys_density * old_parcel_cloud->radius[p_idx] * old_parcel_cloud->radius[p_idx]);

    // tab_omsq is omega^2 in TAB paper referenced in header -- sys_density used in place of rho_l
    CONVERGE_precision_t tab_omsq = spray_tab_csubk * old_parcel_cloud->surf_ten[p_idx] /
                                        (sys_density * old_parcel_cloud->radius[p_idx] *
                                         old_parcel_cloud->radius[p_idx] * old_parcel_cloud->radius[p_idx]) -
                                    tab_rtd * tab_rtd;

    if (tab_omsq <= 0.0)
    {
        // printf("omega squared < 0\n");
        old_parcel_cloud->distort[p_idx] = 0.0;
        old_parcel_cloud->distort_dot[p_idx] = 0.0;
    }
    else
    {
        CONVERGE_precision_t tab_om = CONVERGE_sqrt(tab_omsq);
        CONVERGE_precision_t term1 = old_parcel_cloud->distort[p_idx] - spray_tab_cfocbck * weber;
        CONVERGE_precision_t term2 = (1.0 / tab_om) * (old_parcel_cloud->distort_dot[p_idx] + term1 * tab_rtd);
        // Don't update y and ydot yet, use old values        
        // following is Eq. 4 in TAB paper referenced in header
        CONVERGE_precision_t cos_omdt,sin_omdt;
        cos_omdt = cos(tab_om*dt);
        sin_omdt = sin(tab_om*dt);
        old_parcel_cloud->distort[p_idx] =
            spray_tab_cfocbck * weber + exp(-dt * tab_rtd) * (term1 * cos_omdt + term2 *sin_omdt); // Eq 11.92 in CONVERGE v3.0 Manual
        // printf("y = %f\n",old_parcel_cloud.distort[p_idx]);
        old_parcel_cloud->distort_dot[p_idx] =
            (spray_tab_cfocbck * weber - old_parcel_cloud->distort[p_idx]) * tab_rtd +
            exp(-dt * tab_rtd) * tab_om * (term2 * cos_omdt - term1 * sin_omdt); // Eq 11.93 in CONVERGE
        // make sure y is between 0 and 1
        if (old_parcel_cloud->distort[p_idx] >= 1.0)
        {
            old_parcel_cloud->distort[p_idx] = 1.0;
        }
        if (old_parcel_cloud->distort[p_idx] <= 0.0)
        {
            old_parcel_cloud->distort[p_idx] = 0.0;
        }
    }


}