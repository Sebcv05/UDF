//BreakupCriterion.c
//This function calculates the breakup criterion kb 
#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <BreakupCriterion.h>

CONVERGE_precision_t BreakupCriterion(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx,CONVERGE_precision_t dt)
{


            // Breakup Integral
            // Get value for previous timestep
            CONVERGE_precision_t trap, eta_0, eta, kb,omega,omega_tm1,int_omega;
            omega = old_parcel_cloud->omega[p_idx];
            omega_tm1 = old_parcel_cloud->omega_tm1[p_idx];
            int_omega = old_parcel_cloud->int_omega[p_idx];

            trap = 0.5 * dt * (omega + omega_tm1);      //Trapezium contirbution to integral of omega
            int_omega = int_omega + trap;               //Update Omega

            //printf("\ntrap %e omega_tm1 %e",trap,omega_tm1);
            // eta_0 = 0.05 * old_parcel_cloud->r_drop_0[p_idx];                        // Eta_0 set to a fraction of the initial drop radius
            eta_0 = old_parcel_cloud->eta_drop_0[p_idx];
            // Update
            CONVERGE_precision_t e = 2.71828182845904523536;
            eta = eta_0 * pow(e,int_omega); // Current instability amplitude 
            // Breakup criteria kb
            if(old_parcel_cloud->r_therm[p_idx]<old_parcel_cloud->r_bubble[p_idx])
            {
                //printf("warning, r_bubble > r_therm\n");
            }
   
            kb = eta / (old_parcel_cloud->radius[p_idx] - old_parcel_cloud->r_bubble[p_idx]); // Breakup criteria
            //       if(kb>0.9)
            //    {
            //       //printf("\nt_parcel = %e omega = %e omega_tm1 = %e kb = %e eta_0 = %e eta = %e int_om = %e exp(int_om) = %e",old_parcel_cloud->lifetime[p_idx],old_parcel_cloud->omega[p_idx],old_parcel_cloud->omega_tm1[p_idx],kb,eta_0,eta,int_omega,pow(e,int_omega));
            //    }
           // printf("\nomega = %e omega_tm1 = %e kb = %e eta_0 = %e eta = %e int_om = %e exp(int_om) = %e",old_parcel_cloud->omega[p_idx],old_parcel_cloud->omega_tm1[p_idx],kb,eta_0,eta,int_omega,pow(e,int_omega));
            //Update omega_tm1 & int_omega
            old_parcel_cloud->omega_tm1[p_idx] = omega;
            old_parcel_cloud->int_omega[p_idx] = int_omega;
            // printf("\n kb = %e   eta_0 %e    int_omega %e Rb %e Rd %e  ",kb,eta_0,int_omega,Rb,old_parcel_cloud->radius[p_idx]);
	//      char *filename = "kb_log.txt";
        // FILE *fp = fopen("kb_log.txt","ar");
        // if (fp == NULL)
        //         {
        //         printf("Error opening the file %s",filename);
        //         }
        //  //       fprintf(fp,"\n%i    %e      %f  %e",old_parcel_cloud->identity[p_idx],old_parcel_cloud->lifetime[p_idx],kb,old_parcel_cloud->radius[p_idx]-old_parcel_cloud->r_bubble[p_idx]);
       // fclose(fp);
        
                // printf("\nkb    = %e\n omega = %e   int_omega = %e, eta = %e  t_parcel = %e r_bubble = %e r_drop = %e",kb,old_parcel_cloud->omega[p_idx],int_omega,eta,old_parcel_cloud->lifetime[p_idx],old_parcel_cloud->r_bubble[p_idx],old_parcel_cloud->radius[p_idx]);

return(kb); 
}