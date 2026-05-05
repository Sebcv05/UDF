#ifndef PARCEL_RESET_H
#define PARCEL_RESET_H

#include "lagrangian/env.h"
#include <CONVERGE/udf.h>

/**
 * @brief Reset a parcel to child state with injection radius
 * 
 * This function is called when thermal bubble growth must cease due to:
 * - Subcooling conditions (P_sat < P_amb or T_drop < T_sat)
 * - Bubble collapse (Rdot < 0)
 * - Very small bubble velocity
 * - Droplet radius too small
 * - Other safety/termination conditions
 * 
 * The reset performs:
 * 1. Zeros bubble radius (r_bubble = 0)
 * 2. Resets droplet to injection radius (r_drop_0)
 * 3. Conserves mass by adjusting num_drop: Nd_new = Nd_old * (R_old/R_new)^3
 * 4. Sets pbt = 0 (disables thermal breakup tracking)
 * 5. Sets is_child = 1 (marks as child so KH-RT/evaporation continue)
 * 6. Sets thermal_breakup_flag = 999 (permanent abort flag)
 * 
 * @param parcel_cloud Pointer to ParcelCloud structure
 * @param p_idx Local parcel index
 * @param reason Short string describing why reset occurred (for logging)
 */
void reset_parcel_to_child(struct ParcelCloud* parcel_cloud, 
                            CONVERGE_index_t p_idx, 
                            const char* reason);

#endif // PARCEL_RESET_H
