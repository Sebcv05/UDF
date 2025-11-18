//Vb.c

#include "lagrangian/env.h"
#include <parcel_reset.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <Vb.h>

void Bubble_Velocity(struct ParcelCloud* old_parcel_cloud,CONVERGE_index_t p_idx,CONVERGE_precision_t P_sat, CONVERGE_precision_t P_amb)
{
        CONVERGE_precision_t Vb_temp;
      // Calculate Vb
            if ((P_sat - P_amb) < 0.0)
            {
               // Vb = -pow((3 / rho_l) * (P_amb - P_sat), 0.5);
               //
              // printf("\nVb negative, continuing ...  P_sat = %e P_amb = %e Vb = %e Td = %f t_parcel = %e r_drop = %e\n", P_sat, P_amb, Vb_temp,old_parcel_cloud->temp[p_idx],old_parcel_cloud->lifetime[p_idx],old_parcel_cloud->radius[p_idx]);
               reset_parcel_to_child(old_parcel_cloud, p_idx, "Vb calc: P_sat < P_amb");
               return;
            }
            else
            {
               Vb_temp = CONVERGE_sqrt((3.0 / old_parcel_cloud->density[p_idx]) * (P_sat - P_amb));
               if(Vb_temp<1.0e-10)
               {
                  printf("\nVb temp in bubble velocity function is 0 ");
                  
               }
            }
            
            old_parcel_cloud->v_bubble[p_idx] = Vb_temp ;
}