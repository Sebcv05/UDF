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


/**********************************************************************/
/*                                                                    */
/* Name: user_parcel_prop                                             */
/*                                                                    */
/* Description: user_parcel_prop sets values for the     custom       */
/* parcel properties when new parcels are created. This is            */
/* accomplished through the four functions: user_parcel_inject,       */
/* user_parcel_child, user_parcel_film, user_parcel_splash            */
/*                                                                    */
/*                                                                    */
/**********************************************************************/

// initialize values for the custom parcel properties when new parcels
// are first injected into the domain
// CONVERGE_void_t CONVERGE_register_parcel_inject(CONVERGE_string_t user_name, CONVERGE_ParcelInject_t function)


void custom_parcel_inject(CONVERGE_ParcelInjectInput_t* input)
// CONVERGE_UDF(parcel_inject,
//              IN(VALUE(CONVERGE_index_t, input->parcel_idx), VALUE(CONVERGE_cloud_t, input->spray_cloud)),
//              OUT(CONVERGE_VOID))
{
   // printf("\n START OF PARCEL_PROP.C \n");
   struct ParcelCloud parcel_cloud;

   // printf("PARCEL_PROP.C L58\n");
   load_user_cloud(&parcel_cloud, input->spray_cloud);
   
 
   // printf("PARCEL_PROP.C L70\n");

   // // printf("\n ambient_pres = %f \n",ambient_pres);
   // // printf("\n r_bubble old = %f \n", parcel_cloud.r_bubble[input->parcel_idx]);
   // printf("\n p_idx = %i ", input->parcel_idx);
   // parcel_cloud.radius_bubble[input->parcel_idx] = 1e-10;
   // printf("PARCEL_PROP.C L55\n");



   // // printf("\n PARCEL_PROP.C L69 r_bubble = %e\n", parcel_cloud.r_bubble[input->parcel_idx]);
   // parcel_cloud.v_bubble[input->parcel_idx] = 0.0;
   // // parcel_cloud.r_bubble_0[input->parcel_idx] = parcel_cloud.r_bubble[input->parcel_idx];
   // // printf("\n END OF PARCEL_PROP.C \n");
   // // printf("\n\n r_bubble = %e 	r_bubble_0 = %e", parcel_cloud.r_bubble[input->parcel_idx], parcel_cloud.r_bubble_0[input->parcel_idx]);

   // // R_D_0
   // parcel_cloud.r_drop_0[input->parcel_idx] = parcel_cloud.radius[input->parcel_idx];

   // // // Droplet radius for thermal breakup purposes
   // // parcel_cloud.r_therm[input->parcel_idx] = parcel_cloud.radius[input->parcel_idx];  //this line doesn't work 
   // // // Zero Omega and Eta (breakup variables)
   // parcel_cloud.int_omega[input->parcel_idx] = 0;
   // parcel_cloud.omega[input->parcel_idx] = parcel_cloud.omega_tm1[input->parcel_idx] = 0;
   // parcel_cloud.eta_drop[input->parcel_idx] = 0;
   // parcel_cloud.dgre_cycle_count[input->parcel_idx] = 0;
   // parcel_cloud.tbt[input->parcel_idx] = 0;
   // parcel_cloud.pbt[input->parcel_idx] = 1;
   // parcel_cloud.m0[input->parcel_idx] = (1.33333 * PI * CONVERGE_cube(parcel_cloud.radius[input->parcel_idx]) * parcel_cloud.num_drop[input->parcel_idx]);
   // Set breakup flag to 0
   // parcel_cloud.thermal_breakup_flag[input->parcel_idx] = -1;
   // parcel_cloud.user_parcel_index[input->parcel_idx] = CONVERGE_random_precision();
   parcel_cloud.cloud_index[input->parcel_idx] = CONVERGE_random_precision();
   parcel_cloud.user_parcel_index[input->parcel_idx] = CONVERGE_random_precision();
   // char *filename = "DGRE.txt";
   // char *filename1 = "Temp_Tracker.txt";
   // CONVERGE_index_t ncyc = CONVERGE_ncyc();
   // // Populate Text File With Header
   // if (ncyc == 12)
   // {
   //    FILE *fp = fopen("DGRE.txt", "w");
   //    if (fp == NULL)
   //    {
   //       printf("Error opening the file %s", filename);
   //    }
   //    fprintf(fp, "a   b       c       d       real(x0)      real(x1)      real(x2)      omega   eta	T_drop	P_sat	Rb	Rd     kb	t_parcel	breakup_flag\n");
   //    fclose(fp);
   //    FILE *fp1 = fopen("Temp_Tracker.txt", "w");
   //    if (fp1 == NULL)
   //    {
   //       printf("Error opening the file %s", filename1);
   //    }
   //    fprintf(fp1, "CI    PI    Td    Rd    age    Vmag\n");
   //    fclose(fp1);
   // }
   // printf("\n parcel_prop.c END\n");
}


void custom_parcel_child(CONVERGE_ParcelInjectInput_t* input)
{

}
void custom_parcel_splash(CONVERGE_ParcelSplashInput_t* input)
{

}
void custom_parcel_strip(CONVERGE_ParcelStripInput_t* input)
{

}