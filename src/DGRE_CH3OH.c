//DGRECH3OH.c

#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <DGRE_CH3OH.h>
#include <BubbleDensityCH3OH.h>
#include <PsatCH3OH.h>
#include <CubicSolver.h>

void DGRE_CH3OH(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx,CONVERGE_precision_t global_density)
{
            if(old_parcel_cloud->v_bubble[p_idx]<1.0e-20)
            {
               return;
            }
            // Cycle Count Iterator
            old_parcel_cloud->dgre_cycle_count[p_idx]++;

            //Get Bubble Density Estimate
            CONVERGE_precision_t P_sat;
             Saturation_PressureCH3OH(old_parcel_cloud->temp[p_idx],&P_sat);
            CONVERGE_precision_t rho_b = bubble_densityCH3OH(P_sat,old_parcel_cloud->temp[p_idx]);
            //  Non-dimensionalising terms
            CONVERGE_precision_t Rb, c_omega, We_o, We_i, Ma_i, Delta, psi_o, psi_i,inv_surten,Vdsquared,Vbsquared,inv_Delta, Delta_over_Ma_i,Delta_sq;
            CONVERGE_precision_t sound_speed = 337.44; // approx number for CH3OH at 360 K 2.4 bar (perhaps should consider some function for this as it varies a fair amount with P,T) https://doi.org/10.1016/j.jct.2007.07.002
            Rb = old_parcel_cloud->r_bubble[p_idx];
            Vdsquared = old_parcel_cloud->v_drop[p_idx]*old_parcel_cloud->v_drop[p_idx];
            Vbsquared = old_parcel_cloud->v_bubble[p_idx]*old_parcel_cloud->v_bubble[p_idx];

            inv_surten= 1.0/old_parcel_cloud->surf_ten[p_idx];          //Defined here to reduce total number of division operations 
            c_omega = CONVERGE_sqrt((old_parcel_cloud->density[p_idx] * Rb*Rb*Rb ) * inv_surten);
            We_o = old_parcel_cloud->density[p_idx] * Vdsquared * Rb * inv_surten;
            We_i = old_parcel_cloud->density[p_idx] * Vbsquared * Rb * inv_surten;
            Ma_i = old_parcel_cloud->v_bubble[p_idx] / sound_speed;
            Delta = old_parcel_cloud->radius[p_idx] / old_parcel_cloud->r_bubble[p_idx];
            Delta_sq = Delta * Delta;
            inv_Delta = 1.0/Delta;
            psi_o = global_density / old_parcel_cloud->density[p_idx];
            psi_i = rho_b / old_parcel_cloud->density[p_idx];
            //printf("\nWe_o = %e We_i = %e Ma_i = %e Delta = %e psi_o = %e psi_i = %e, n = %i", We_o, We_i, Ma_i, Delta, psi_o, psi_i, old_parcel_cloud->dgre_cycle_count[p_idx]);
            // See documentation on DGRE for definitions of alpha-epsilon
            CONVERGE_precision_t dgre_alpha, dgre_beta, dgre_gamma, dgre_delta, dgre_epsilon;
            dgre_alpha = Delta - Delta_sq - psi_o * Delta;
            dgre_beta = CONVERGE_sqrt(We_o) * ((Delta_sq*Delta_sq) + psi_o - 1.0);
            dgre_gamma = 2.0 * (inv_Delta*inv_Delta + Delta_sq);
            Delta_over_Ma_i = Delta/Ma_i;
            dgre_delta = 3.0 * psi_i * We_i * Delta_over_Ma_i * Delta_over_Ma_i;
            dgre_epsilon = 3.0 * CONVERGE_sqrt(We_i);
            // See latex documentation on DGRE
            CONVERGE_precision_t a, b, c, d;
            a = dgre_alpha;
            b = dgre_beta + dgre_alpha * dgre_epsilon;
            c = dgre_beta * dgre_epsilon + dgre_gamma - dgre_delta;
            d = dgre_gamma * dgre_epsilon;
            //printf("\na = %e b =  %e c =  %e d = %e ", a, b, c, d);
            // WARNING FOR NOT A NUMBER ERRROR
            if (isnan(a) || isnan(b) || isnan(c) || isnan(d))
            {
               printf("ERROR - DGRE CO-EFFICIENTS ARE NAN");
               printf("\nWe_i	%e	We_o 	%e	Delta	%e	Psi_i	%e	Psi_o	%e	c_omega	%e n %i\n", We_i, We_o, Delta, psi_i, psi_o, c_omega, old_parcel_cloud->dgre_cycle_count[p_idx]);
               printf("continuing on L471");
               return;
            }

            struct roots r1, r2;

            r1 = Cubic_Solver(a, b, c, d);
            // Unit test for Cubic Solver

            // Manual Unit Testing --Disabled for speed at runtime 
           /* r2 = Cubic_Solver(5, 3, 2, 8);
            assert(creal(r2.x0) == -1.27305); /// Don't know why but this does not do anyting
            if (roundf(creal(r2.x0) * 1e5) / 1e5 != -1.27305)
            {
               printf("\nCubic Solver has Failed Unit Test!");
               printf("-1.27305 != %f", creal(r2.x0));
               CONVERGE_mpi_abort();
            }
            r2 = Cubic_Solver(5, 8, 2, -8);
            if (roundf(creal(r2.x0) * 1e5) / 1e5 != -1.17255)
            {
               printf("\nCubic Solver has Failed Unit Test!");
               printf("-1.17255 != %f", creal(r2.x0));
               CONVERGE_mpi_abort();
            }
            r2 = Cubic_Solver(0, 8, 2, 8);
            if (roundf(creal(r2.x0) * 1e5) / 1e5 != -0.125)
            {
               printf("\nCubic Solver has Failed Unit Test!");
               printf("-0.125 != %f", creal(r2.x0));
               CONVERGE_mpi_abort();
            }
            */
           
            // Select root with largest real part (maximum disturbance growth rate)
            CONVERGE_precision_t root_values[3];
            root_values[0] = creal(r1.x0);
            root_values[1] = creal(r1.x1);
            root_values[2] = creal(r1.x2);
            
            // Diagnostic: Print all roots for first few calls
            static int root_diag_count = 0;
            if (root_diag_count < 5) {
               printf("[ROOT_DIAG] roots: x0=%.3e, x1=%.3e, x2=%.3e\n", 
                      root_values[0], root_values[1], root_values[2]);
               root_diag_count++;
            }
            
            // Find root with maximum real part
            CONVERGE_precision_t max_real = root_values[0];
            if (root_values[1] > max_real) max_real = root_values[1];
            if (root_values[2] > max_real) max_real = root_values[2];
            
            // Check for unreasonably large values
            if (fabs(max_real) > 1.0e10)
            {
               max_real = 0.0;
            }
            
            // If negative, use absolute value (modulus)
            if (max_real < 0.0)
            {
               max_real = fabs(max_real);
            }
            
            // Convert from non-dimensional to dimensional omega
            old_parcel_cloud->omega[p_idx] = max_real / c_omega;
}