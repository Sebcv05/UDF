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

/** Load property tables for all phase species in iterator
 */
void load_species_tables(CONVERGE_iterator_t it, CONVERGE_index_t table_id, CONVERGE_table_t **tables)
{
   // Allocate the number of tables denoted my the species iterator
   *tables = (CONVERGE_table_t *)malloc(CONVERGE_iterator_distance(it) * sizeof(CONVERGE_table_t));

   for(CONVERGE_index_t local_isp = 0, i = CONVERGE_iterator_first(it); i != -1;
       local_isp++, i                    = CONVERGE_iterator_next(it))
   {
      (*tables)[local_isp] = CONVERGE_get_species_prop_table(table_id, i);
   }
}

/** Unload property tables for all phase species in iterator
 */
void unload_species_tables(CONVERGE_iterator_t it, CONVERGE_table_t **tables)
{
   // Cleanup local references to internal tables
   for(CONVERGE_index_t local_isp = 0, i = CONVERGE_iterator_first(it); i != -1;
       local_isp++, i                    = CONVERGE_iterator_next(it))
   {
      CONVERGE_table_destroy((*tables) + local_isp);
   }
   free(*tables);
}

/** Compute the parcel property using a table
 */
double get_parcel_prop_from_table(CONVERGE_iterator_t it,
                                  CONVERGE_table_t **table,
                                  CONVERGE_precision_t temp,
                                  CONVERGE_precision_t *mfrac,
                                  CONVERGE_index_t is_density)
{
   double prop_value = 0.0;

   // If single component, no need to do mass-averaging.
   if(CONVERGE_iterator_distance(it) == 1)
   {
      prop_value = CONVERGE_table_lookup(*table[0], temp);
      return prop_value;
   }

   for(CONVERGE_index_t local_isp = 0, i = CONVERGE_iterator_first(it); i != -1;
       local_isp++, i                    = CONVERGE_iterator_next(it))
   {
      if(is_density)
      {
         prop_value += mfrac[local_isp] / CONVERGE_table_lookup(*table[local_isp], temp);
      }
      else
      {
         prop_value += mfrac[local_isp] * CONVERGE_table_lookup(*table[local_isp], temp);
      }
   }

   if(is_density)
   {
      prop_value = 1.0 / prop_value;
   }

   return prop_value;
}
