// Geometry.c

#include "lagrangian/env.h"
#include <user_header.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <Geometry.h>

void Geometry(struct ParcelCloud *old_parcel_cloud, CONVERGE_index_t p_idx, CONVERGE_precision_t dt)
{
   // UPDATE DROPLET PROPERTIES USING GEOMETRIC ARGUMENT DELTA^2 = Vd/Vb
   CONVERGE_precision_t r_therm_tm1=old_parcel_cloud->r_therm[p_idx];
   CONVERGE_precision_t old_rad = old_parcel_cloud->radius[p_idx];
   CONVERGE_precision_t old_nd = old_parcel_cloud->num_drop[p_idx];
   CONVERGE_precision_t Vd = old_parcel_cloud->v_bubble[p_idx] * pow(old_parcel_cloud->radius[p_idx] / old_parcel_cloud->r_bubble[p_idx], -2);
   if (Vd < 0.0)
   {
      printf("\n Droplet imploding, aborting\n");
      CONVERGE_mpi_abort();
   }
   old_parcel_cloud->v_drop[p_idx] = Vd;
   CONVERGE_precision_t dRd = Vd * dt;
   if (dRd < 0.0)
   {
      printf("\ndRd<0,setting tbt to 999 and continuing");
      old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
      return;
   }
   if (old_parcel_cloud->r_therm[p_idx] + dRd > 0 && dRd < old_parcel_cloud->r_therm[p_idx])
   {
      // printf("\nupdating r_therm (less than double)");
      //  old_parcel_cloud->r_therm[p_idx]
      CONVERGE_precision_t r_therm_temp = old_parcel_cloud->radius[p_idx] + dRd;
      if (r_therm_temp < 1e-10)
      {
         printf("\n r_therm being set to 0 r_therm prev = %e dRd = %e", old_parcel_cloud->r_therm[p_idx], dRd);
         CONVERGE_mpi_abort();
      }
      else if (r_therm_temp >= old_parcel_cloud->radius[p_idx])
      {
         // old_parcel_cloud->r_therm[p_idx] = r_therm_temp;
         old_parcel_cloud->radius[p_idx] = r_therm_temp;
         old_parcel_cloud->num_drop[p_idx] = old_nd * CONVERGE_cube(old_rad/old_parcel_cloud->radius[p_idx]);
      }
      else // 0<r_therm <r_drop         -- This case shouldn't happen as dRd>0
      {
         printf("\n Geometry.c");
         printf("\nr_therm smaller than r_drop,aborting");
         printf("\n r_therm = %e r_therm_tm1 = %e r_drop = %e, rb= %e vb= %e vd = %e dRd = %e",old_parcel_cloud->r_therm[p_idx],r_therm_tm1,old_parcel_cloud->radius[p_idx],old_parcel_cloud->r_bubble[p_idx],old_parcel_cloud->v_bubble[p_idx],Vd,dRd);
         CONVERGE_mpi_abort();
      }
      // printf("\n r_therm = %e radius = %e",old_parcel_cloud->r_therm[p_idx],old_parcel_cloud->radius[p_idx]);
   }
   else if (dRd > old_parcel_cloud->radius[p_idx]) // This if statement has no physical justification other than the assumption that doubling the droplet's radius in dt will cause breakup
   {
      // printf("\nGeometry.c\n");
      // //printf("\ntriggering breakup because geometry leads to doubling of droplet radius");
      // printf("\nAborting because Geometry.c has doubled droplet radius in dt.....");
      // printf("\n Vb = %e Rb = %e Rd_old = %e Vd = %e Rb_new = %e",old_parcel_cloud->v_bubble[p_idx],old_parcel_cloud->r_bubble[p_idx],old_parcel_cloud->r_therm[p_idx],Vd,old_parcel_cloud->r_therm[p_idx]+dRd);

      // CONVERGE_mpi_abort();
      old_parcel_cloud->tbt[p_idx] =1;
      return;
   }
   if (old_parcel_cloud->radius[p_idx] < 1.0e-10)
   {
      printf("\nvery small radius = %e\n", old_parcel_cloud->r_therm[p_idx]);
   }
}