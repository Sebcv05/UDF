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
   
   // SAFETY: Prevent division by very small r_bubble
   CONVERGE_precision_t r_bubble_safe = old_parcel_cloud->r_bubble[p_idx];
   if (r_bubble_safe < 1e-12) {
      r_bubble_safe = 1e-12;
   }
   
   CONVERGE_precision_t Vd = old_parcel_cloud->v_bubble[p_idx] * pow(old_parcel_cloud->radius[p_idx] / r_bubble_safe, -2);
   
   // SAFETY: Prevent negative Vd (droplet shrinkage during thermal breakup)
   if (Vd < 0.0)
   {
      static int negative_vd_count = 0;
      if (negative_vd_count < 3) {
         printf("\n[GEOMETRY] Negative Vd=%.3e, v_bubble=%.3e, setting Vd=0 and continuing\n", 
                Vd, old_parcel_cloud->v_bubble[p_idx]);
         negative_vd_count++;
      }
      Vd = 0.0;  // Don't shrink, just keep current size
   }
   
   old_parcel_cloud->v_drop[p_idx] = Vd;
   CONVERGE_precision_t dRd = Vd * dt;
   
   if (dRd < 0.0)
   {
      printf("\n[GEOMETRY] dRd<0, setting tbt to 999 and continuing");
      old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
      return;
   }
   if (old_parcel_cloud->r_therm[p_idx] + dRd > 0 && dRd < old_parcel_cloud->radius[p_idx])
   {
      // printf("\nupdating r_therm (less than double)");
      //  old_parcel_cloud->r_therm[p_idx]
      CONVERGE_precision_t r_therm_temp = old_parcel_cloud->radius[p_idx] + dRd;
      if (r_therm_temp < 1e-10)
      {
         printf("\n[GEOMETRY] r_therm would be too small (%.3e), keeping current radius\n", r_therm_temp);
         // Don't update, keep current radius
         return;
      }
      else if (r_therm_temp >= old_parcel_cloud->radius[p_idx])
      {
         // Droplet growing - update radius and num_drop
         // SAFETY: Don't allow radius to shrink below original
         CONVERGE_precision_t new_radius = r_therm_temp;
         if (new_radius < old_rad) {
            new_radius = old_rad;  // Keep at least original size
         }
         
         old_parcel_cloud->radius[p_idx] = new_radius;
         old_parcel_cloud->num_drop[p_idx] = old_nd * CONVERGE_cube(old_rad/new_radius);
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