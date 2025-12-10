//RPE_euler.c
// Rayleigh-Plesset Equation solver with Explicit Euler integration

/*
 * breakup_phase states:
 *   0 = DISABLED  (parent, not eligible - subcooled, too small, etc.)
 *   1 = ELIGIBLE  (parent, superheated, ready to enter thermal breakup)
 *   2 = ACTIVE    (parent, growing bubble in sub-timestep loop)
 *   3 = RECOVERY  (parent, bubble collapsed, attempting recovery)
 *   4 = READY     (parent, bubble at threshold, ready to fragment)
 *   5 = COMPLETE  (child - result of actual breakup)
 *   
 *   DIAGNOSTIC STATES (bypassed breakup):
 *   12 = Droplet too small (RPE_euler line 258)
 *   13 = P_sat < P_amb subcooled (RPE_euler line 339)
 *   14 = Recovered parcel in RPE (RPE_euler line 375)
 *   15 = Bubble collapse Rdot < 0 (RPE_euler line 393)
 *   16 = Subcooled T < T_sat (RPE_euler line 426)
 *   17 = Rdot too small (RPE_euler line 460)
 */

#include "lagrangian/env.h"
#include <RPE_euler.h>
#include <PsatNH3.h>
#include <TsatNH3.h>
#include <BubbleDensityNH3.h>
#include <parcel_reset.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>

// Physical constants
#define PI 3.14159265358979323846
#define R_SPEC_NH3 488.2  // J/(kg·K) - Specific gas constant for ammonia

// Debug logging configuration
#define RPE_DEBUG_LOGGING 1        // Set to 0 to disable
#define LOG_INTERVAL 50            // Log every Nth call
#define LOG_FILE_NAME "rpe_parcel_debug.csv"

// Utility functions
CONVERGE_precision_t safe_sqrt(CONVERGE_precision_t x) {
    return (x > 0.0) ? CONVERGE_sqrt(x) : 0.0;
}

CONVERGE_precision_t safe_divide(CONVERGE_precision_t num, CONVERGE_precision_t denom, CONVERGE_precision_t fallback) {
    return (fabs(denom) > 1e-20) ? (num / denom) : fallback;
}

// Compute heat and mass transfer rates
void compute_thermal_mass_transfer(
    CONVERGE_precision_t R,
    CONVERGE_precision_t Rdot,
    CONVERGE_precision_t T_drop,
    CONVERGE_precision_t m_b,
    const RPE_Params* params,
    CONVERGE_precision_t* Nu_out,
    CONVERGE_precision_t* Q_out,
    CONVERGE_precision_t* mdot_out
) {
    // Calculate bubble pressure from current mass
    CONVERGE_precision_t Vb = (4.0/3.0) * PI * R*R*R;
    CONVERGE_precision_t rho_v = safe_divide(m_b, Vb, 1e-6);
    if (rho_v < 1e-6) rho_v = 1e-6;
    
    CONVERGE_precision_t T_drop_safe = (T_drop > 1e-3) ? T_drop : 1e-3;
    CONVERGE_precision_t Pb = rho_v * params->R_spec * T_drop_safe;
    
    // Interface temperature from bubble pressure
    CONVERGE_precision_t T_int = T_satNH3(Pb);
    
    // Temperature difference for heat transfer
    CONVERGE_precision_t dT = T_drop - T_int;
    
    // Reynolds and Prandtl numbers
    CONVERGE_precision_t L_char = 2.0 * R;
    CONVERGE_precision_t Re = params->rho_l * fabs(Rdot) * L_char / params->mu_l;
    CONVERGE_precision_t Pr = params->mu_l * params->cp_l / params->k_l;
    
    // Nusselt number (Ranz-Marshall correlation)
    CONVERGE_precision_t Nu_uncapped = 2.0 + 0.6 * safe_sqrt(Re) * pow(Pr, 1.0/3.0);
    CONVERGE_precision_t Nu = (Nu_uncapped < params->max_Nu) ? Nu_uncapped : params->max_Nu;
    
    // Surface area
    CONVERGE_precision_t A_bubble = 4.0 * PI * R * R;
    
    // ========== DUAL GEOMETRY: Semi-infinite vs Film ==========
    // Match Python euler_explicit.py: evaluate both, use max
    
    // Semi-infinite formulation
    CONVERGE_precision_t h_conv_si = Nu * params->k_l / L_char;
    CONVERGE_precision_t Q_si = h_conv_si * A_bubble * dT;
    
    // Film formulation
    CONVERGE_precision_t film_thickness = params->Ro - R;
    CONVERGE_precision_t h_conv_film, Q_film;
    
    if (film_thickness > 1e-12) {
        h_conv_film = Nu * params->k_l / film_thickness;
        Q_film = h_conv_film * A_bubble * dT;
    } else {
        h_conv_film = 1e20;
        Q_film = 0.0;
    }
    
    // Use maximum (least restrictive) for heat transfer
    CONVERGE_precision_t Q_conv = (Q_film > Q_si) ? Q_film : Q_si;
    
    // Mass transfer rate (thermal limiting)
    CONVERGE_precision_t mdot = safe_divide(Q_conv, params->L_v, 0.0);
    if (mdot < 0.0) mdot = 0.0;  // Evaporation only
    
    // Output
    *Nu_out = Nu;
    *Q_out = Q_conv;
    *mdot_out = mdot;
}

// Compute derivatives for the 4 ODEs
void compute_derivatives(
    const BubbleState* state,
    const RPE_Params* params,
    BubbleDerivatives* derivs,
    CONVERGE_precision_t mdot
) {
    CONVERGE_precision_t R = state->R;
    CONVERGE_precision_t Rdot = state->Rdot;
    CONVERGE_precision_t T_drop = state->T_drop;
    CONVERGE_precision_t m_b = state->m_b;
    
    // Safety checks
    if (R < 1e-12) R = 1e-12;
    
    // Calculate bubble pressure
    CONVERGE_precision_t Vb = (4.0/3.0) * PI * R*R*R;
    CONVERGE_precision_t rho_v = m_b / Vb;  // Direct calculation - Vb cannot be zero if R > 0
    if (rho_v < 1e-6) rho_v = 1e-6;  // Minimum density floor
    
    CONVERGE_precision_t T_drop_safe = (T_drop > 1e-3) ? T_drop : 1e-3;
    CONVERGE_precision_t Pb = rho_v * params->R_spec * T_drop_safe;
    
    // DIAGNOSTIC: Check bubble pressure calculation (first few calls)
    static int pb_diagnostic_count = 0;
    if (pb_diagnostic_count < 3 && R < 1e-7) {
        printf("[RPE_PB_DIAG] R=%.3e m, m_b=%.3e kg, Vb=%.3e m³, rho_v=%.3e kg/m³\n",
               R, m_b, Vb, rho_v);
        printf("              Pb=%.3e Pa, P_amb=%.3e Pa, DeltaP=%.3e Pa\n",
               Pb, params->P_amb, Pb - params->P_amb);
        pb_diagnostic_count++;
    }
    
    // 1. dR/dt = Rdot
    derivs->dRdt = Rdot;
    
    // 2. Rayleigh-Plesset equation: dRdot/dt
    CONVERGE_precision_t pressure_term = safe_divide(Pb - params->P_amb - 2.0*params->sigma/R - 4.0*params->mu_l*Rdot/R,
                                                       params->rho_l * R, 0.0);
    CONVERGE_precision_t inertial_term = -1.5 * Rdot * Rdot / R;
    CONVERGE_precision_t Rddot = pressure_term + inertial_term;
    
    // Apply acceleration limits
    if (Rddot > 1e10) Rddot = 1e10;
    if (Rddot < -1e10) Rddot = -1e10;
    
    derivs->dRdotdt = Rddot;
    
    // Compute current mass transfer rate
    CONVERGE_precision_t Nu, Q_conv, mdot_current;
    compute_thermal_mass_transfer(R, Rdot, T_drop, m_b, params, &Nu, &Q_conv, &mdot_current);
    
    // 3. Mass balance: dm_b/dt = mdot (use current value)
    derivs->dmbdt = mdot_current;
    
    // 4. Energy balance: dT_drop/dt (use same mdot for consistency)
    CONVERGE_precision_t Q_evap = params->L_v * mdot_current;
    CONVERGE_precision_t dTdt = safe_divide(Q_conv - Q_evap, params->m_drop * params->cp_l, 0.0);
    
    // DIAGNOSTIC: Check energy balance (first few calls for parent parcels)
    static int energy_diag_count = 0;
    if (energy_diag_count < 5 && R > 1e-6 && R < 1e-4) {
        printf("[RPE_ENERGY] R=%.3e m, T_drop=%.2f K, m_drop=%.3e kg\n", R, T_drop, params->m_drop);
        printf("             Nu=%.2f, Q_conv=%.3e W, mdot=%.3e kg/s, Q_evap=%.3e W\n", 
               Nu, Q_conv, mdot_current, Q_evap);
        printf("             Q_net=%.3e W, dTdt=%.3e K/s\n", Q_conv - Q_evap, dTdt);
        energy_diag_count++;
    }
    
    // Apply temperature rate limits
    if (dTdt > 1e6) dTdt = 1e6;
    if (dTdt < -1e6) dTdt = -1e6;
    
    derivs->dTdt = dTdt;
}

// Explicit Euler integration step
void euler_step(
    BubbleState* state,
    const BubbleDerivatives* derivs,
    CONVERGE_precision_t dt
) {
    state->R += derivs->dRdt * dt;
    state->Rdot += derivs->dRdotdt * dt;
    state->T_drop += derivs->dTdt * dt;
    state->m_b += derivs->dmbdt * dt;
    
    // Apply safety limits
    if (state->R < 1e-12) state->R = 1e-12;
    if (state->m_b < 0.0) state->m_b = 0.0;
    if (state->T_drop < 100.0) state->T_drop = 100.0;
    if (state->T_drop > 500.0) state->T_drop = 500.0;
}

// Main solver function - replaces Bubble_Velocity()
void RPE_euler_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
) {
    // One-time logging to confirm thermal model is active
    static int thermal_model_logged = 0;
    if (!thermal_model_logged) {
        printf("\n========================================\n");
        printf("[RPE_MODEL_THERMAL] Thermal RPE model active\n");
        printf("========================================\n\n");
        thermal_model_logged = 1;
    }
    
    // Check if parcel is in recovery mode and print status
    if (old_parcel_cloud->recovery_time[p_idx] > 0.0) {
        CONVERGE_precision_t current_time = CONVERGE_simulation_time_sec();
        CONVERGE_precision_t time_since_recovery = current_time - old_parcel_cloud->recovery_time[p_idx];
        static int recovery_status_count = 0;
        if (recovery_status_count < 20) {
            printf("[RPE_IN_RECOVERY] p_idx=%li, recovery_count=%d, time_since_recovery=%.3e s, R=%.3e m\n",
                   p_idx, old_parcel_cloud->recovery_count[p_idx], time_since_recovery, 
                   old_parcel_cloud->r_bubble[p_idx]);
            recovery_status_count++;
        }
    }
    
    // recovery_time stores the time when last recovery occurred (seconds)
    // recovery_count stores the number of recovery attempts
    
    // Initialize parameters structure
    RPE_Params params;
    
    // Read liquid properties from parcel
    params.rho_l = old_parcel_cloud->density[p_idx];
    params.mu_l = old_parcel_cloud->viscosity[p_idx];
    params.sigma = old_parcel_cloud->surf_ten[p_idx];
    params.P_amb = P_amb;
    params.Ro = old_parcel_cloud->radius[p_idx];
    
    // Safety check: If droplet radius is too small, don't run RPE
    if (params.Ro < 1.0e-9) {
        if (old_parcel_cloud->recovery_time[p_idx] > 0.0) {
            printf("[RPE_KILL_IN_RECOVERY] p_idx=%li, Reason: Ro too small (%.3e m), recovery_time=%.3e s, recovery_count=%d\n",
                   p_idx, params.Ro, old_parcel_cloud->recovery_time[p_idx], old_parcel_cloud->recovery_count[p_idx]);
        }
        printf("[RPE_ERROR] Droplet radius too small: Ro=%.3e m\n", params.Ro);
        old_parcel_cloud->breakup_phase[p_idx] = 12;  // Droplet too small
        old_parcel_cloud->film_flag[p_idx] = 12;
        old_parcel_cloud->r_drop_0[p_idx] = old_parcel_cloud->radius[p_idx];
        old_parcel_cloud->r_bubble[p_idx] = 0.0;
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        return;
    }
    
    // Calculate droplet mass
    params.m_drop = (4.0/3.0) * PI * params.Ro * params.Ro * params.Ro * params.rho_l;
    
    // Get temperature-dependent properties from tables
    CONVERGE_precision_t Td = old_parcel_cloud->temp[p_idx];
    
    // Special handling for parcels re-entering after recovery period
    // These have been reset to breakup_phase=1 and r_bubble=0 by spray_drop_distort
    if (old_parcel_cloud->recovery_time[p_idx] > 0.0 && 
        old_parcel_cloud->r_bubble[p_idx] < 1e-12) {
        
        // Check if conditions are favorable for re-growth
        CONVERGE_precision_t P_sat_check;
        Saturation_PressureNH3(Td, &P_sat_check);
        
        if (P_sat_check > P_amb) {
            // Conditions favorable - initialize new bubble from critical radius
            static int recovery_restart_count = 0;
            if (recovery_restart_count < 10) {
                printf("[RECOVERY_RESTART] p_idx=%li, P_sat=%.3e > P_amb=%.3e, initializing new bubble\n",
                       p_idx, P_sat_check, P_amb);
                printf("                   T_drop=%.2f K, recovery_count=%d\n",
                       Td, old_parcel_cloud->recovery_count[p_idx]);
                recovery_restart_count++;
            }
            
            // Calculate critical radius: Rc = 2*sigma / (P_sat - P_amb)
            CONVERGE_precision_t sigma = old_parcel_cloud->surf_ten[p_idx];
            CONVERGE_precision_t delta_P = P_sat_check - P_amb;
            CONVERGE_precision_t Rc = 2.0 * sigma / delta_P;
            
            // Initialize bubble at 1.1 * Rc (10% above critical for stable growth)
            CONVERGE_precision_t R_init = 1.1 * Rc;
            old_parcel_cloud->r_bubble[p_idx] = R_init;
            old_parcel_cloud->r_bubble_0[p_idx] = R_init;
            old_parcel_cloud->v_bubble[p_idx] = 0.0;
            
            // Initialize bubble mass from saturation properties
            CONVERGE_precision_t rho_b = bubble_densityNH3(P_sat_check, Td);
            CONVERGE_precision_t Vb = (4.0/3.0) * PI * R_init * R_init * R_init;
            old_parcel_cloud->m_bubble[p_idx] = rho_b * Vb;
            
            // Clear recovery time to allow normal operation
            old_parcel_cloud->recovery_time[p_idx] = 0.0;
            
            if (recovery_restart_count < 10) {
                printf("                   Rc=%.3e m, R_init=%.3e m (1.1*Rc), m_bubble=%.3e kg\n",
                       Rc, R_init, old_parcel_cloud->m_bubble[p_idx]);
            }
            
        } else {
            // Conditions not favorable - disable thermal breakup
            static int recovery_abort_count = 0;
            if (recovery_abort_count < 10) {
                printf("[RECOVERY_ABORT] p_idx=%li, P_sat=%.3e < P_amb=%.3e, aborting thermal breakup\n",
                       p_idx, P_sat_check, P_amb);
                recovery_abort_count++;
            }
            old_parcel_cloud->breakup_phase[p_idx] = 13;  // Subcooled
            old_parcel_cloud->film_flag[p_idx] = 13;
            old_parcel_cloud->r_drop_0[p_idx] = old_parcel_cloud->radius[p_idx];
            old_parcel_cloud->r_bubble[p_idx] = 0.0;
            old_parcel_cloud->v_bubble[p_idx] = 0.0;
            return;
        }
    }
    
    // Latent heat and specific heat (weighted by species mass fraction)
    CONVERGE_precision_t L_v_avg = 0.0;
    CONVERGE_precision_t cp_l_avg = 0.0;
    for (CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++) {
        L_v_avg += old_parcel_cloud->mfrac[p_idx * num_parcel_species + isp] * 
                   CONVERGE_table_lookup(hvap_table[isp], Td);
        cp_l_avg += old_parcel_cloud->mfrac[p_idx * num_parcel_species + isp] * 
                    CONVERGE_table_lookup(cp_table[isp], Td);
    }
    
    params.L_v = L_v_avg;
    params.cp_l = cp_l_avg;
    
    // Ammonia-specific properties (TODO: make these configurable)
    params.R_spec = R_SPEC_NH3;
    params.k_l = 0.5;  // Approximate thermal conductivity W/(m·K)
    params.max_Nu = 1000.0;  // Maximum Nusselt number
    
    // Initialize bubble state from parcel
    BubbleState state;
    state.R = old_parcel_cloud->r_bubble[p_idx];
    state.Rdot = old_parcel_cloud->v_bubble[p_idx];
    state.T_drop = Td;
    
    // Initialize or read bubble mass from parcel
    // If m_bubble is zero or very small, initialize from equilibrium
    if (old_parcel_cloud->m_bubble[p_idx] < 1e-20) {
        CONVERGE_precision_t P_sat;
        Saturation_PressureNH3(Td, &P_sat);
        CONVERGE_precision_t rho_b = bubble_densityNH3(P_sat, Td);
        CONVERGE_precision_t Vb = (4.0/3.0) * PI * state.R * state.R * state.R;
        state.m_b = rho_b * Vb;
        old_parcel_cloud->m_bubble[p_idx] = state.m_b;  // Store initial value
    } else {
        // Use stored bubble mass (preserves thermal limiting effect)
        state.m_b = old_parcel_cloud->m_bubble[p_idx];
    }
    
    // Diagnostic: Track RPE solver calls for parent parcels
    static int rpe_call_count = 0;
    static int last_diagnostic_cycle = -1;
    static CONVERGE_precision_t max_R_seen = 0.0;
    static CONVERGE_precision_t max_Rdot_seen = 0.0;
    int current_cycle = CONVERGE_ncyc();
    CONVERGE_precision_t sim_time = CONVERGE_simulation_time();
    
    rpe_call_count++;
    if (state.R > max_R_seen) max_R_seen = state.R;
    if (state.Rdot > max_Rdot_seen) max_Rdot_seen = state.Rdot;
    
    // Report every 100 cycles
    if (current_cycle > 0 && current_cycle % 100 == 0 && current_cycle != last_diagnostic_cycle) {
        CONVERGE_precision_t T_sat_at_P_amb = T_satNH3(P_amb);
        CONVERGE_precision_t P_sat_at_Td;
        Saturation_PressureNH3(Td, &P_sat_at_Td);
        printf("[RPE_STATUS] Cycle %d, Time %.6e s: %d RPE calls, Max R=%.3e µm, Max Rdot=%.3f m/s, T_drop=%.2f K, T_sat=%.2f K, P_sat=%.2e Pa\n",
               current_cycle, sim_time, rpe_call_count, max_R_seen*1e6, max_Rdot_seen, Td, T_sat_at_P_amb, P_sat_at_Td);
        last_diagnostic_cycle = current_cycle;
        rpe_call_count = 0;
        max_R_seen = 0.0;
        max_Rdot_seen = 0.0;
    }
    
    // Check for negative pressure difference
    CONVERGE_precision_t P_sat;
    Saturation_PressureNH3(Td, &P_sat);
    if ((P_sat - P_amb) < 0.0) {
        if (old_parcel_cloud->recovery_time[p_idx] > 0.0) {
            printf("[RPE_KILL_IN_RECOVERY] p_idx=%li, Reason: Negative P_sat-P_amb (%.3e Pa), recovery_time=%.3e s, recovery_count=%d\n",
                   p_idx, (P_sat - P_amb), old_parcel_cloud->recovery_time[p_idx], old_parcel_cloud->recovery_count[p_idx]);
        }
        old_parcel_cloud->breakup_phase[p_idx] = 13;  // P_sat < P_amb subcooled
        old_parcel_cloud->film_flag[p_idx] = 13;
        old_parcel_cloud->r_drop_0[p_idx] = old_parcel_cloud->radius[p_idx];
        old_parcel_cloud->r_bubble[p_idx] = 0.0;
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        return;
    }
    
    // Compute derivatives
    BubbleDerivatives derivs;
    CONVERGE_precision_t Nu, Q_conv, mdot;
    compute_thermal_mass_transfer(state.R, state.Rdot, state.T_drop, state.m_b, 
                                   &params, &Nu, &Q_conv, &mdot);
    compute_derivatives(&state, &params, &derivs, mdot);
    
    // Store before state for diagnostics
    CONVERGE_precision_t R_before = state.R;
    CONVERGE_precision_t Rdot_before = state.Rdot;
    
    // Euler step
    euler_step(&state, &derivs, dt_sub);
    
    // SAFETY: Prevent negative Rdot (bubble collapse) with TIME-BASED recovery logic
    if (state.Rdot < 0.0) {
        // Calculate saturation pressure
        CONVERGE_precision_t P_sat_calc;
        Saturation_PressureNH3(state.T_drop, &P_sat_calc);
        CONVERGE_precision_t T_sat_calc = T_satNH3(params.P_amb);
        
        // Get current simulation time
        CONVERGE_precision_t current_time = CONVERGE_simulation_time_sec();
        
        // Check if this parcel has already been recovered (now a child)
        CONVERGE_precision_t recovery_start_time = old_parcel_cloud->recovery_time[p_idx];
        
        // if (recovery_start_time > 0.0) {
        //     // This parcel was already recovered and converted to child (breakup_phase=5)
        //     // It shouldn't be in RPE anymore
        //     printf("[RPE_ERROR] Recovered parcel (child) re-entered RPE! p_idx=%li, recovery_time=%.3e s\n",
        //            p_idx, recovery_start_time);
        //     old_parcel_cloud->breakup_phase[p_idx] = 14;  // Recovered parcel in RPE
        //     old_parcel_cloud->film_flag[p_idx] = 14;
        //     old_parcel_cloud->r_drop_0[p_idx] = old_parcel_cloud->radius[p_idx];
        //     old_parcel_cloud->r_bubble[p_idx] = 0.0;
        //     old_parcel_cloud->v_bubble[p_idx] = 0.0;
        //     return;
        // }
        
        // First time collapse - attempt recovery by converting to child
        // Always print collapse diagnostics
        printf("[RPE_COLLAPSE] Negative Rdot=%.3e, entering recovery...\n", state.Rdot);
        printf("               Time: %.6e s\n", current_time);
        
        printf("               T_drop=%.2f K, T_sat(P_amb)=%.2f K, P_sat(T_drop)=%.3e Pa, P_amb=%.3e Pa\n",
               state.T_drop, T_sat_calc, P_sat_calc, params.P_amb);
        printf("               R=%.3e m, Ro=%.3e m, dRdt=%.3e m/s, dRdotdt=%.3e m/s²\n",
               state.R, params.Ro, derivs.dRdt, derivs.dRdotdt);
        
        // NEW RECOVERY STRATEGY: Set to RECOVERY state (3) with timed recovery period
        // After recovery period, parcel will be re-evaluated for thermal breakup eligibility
        
        // Apply recovery: reset bubble and enter recovery state
        old_parcel_cloud->breakup_phase[p_idx] = 3;  // RECOVERY (bubble collapsed)
        old_parcel_cloud->film_flag[p_idx] = 3;
        //Use volume balance to return to original N0 
        old_parcel_cloud->num_drop[p_idx] = CONVERGE_cube(old_parcel_cloud->radius[p_idx])* old_parcel_cloud->num_drop[p_idx] / CONVERGE_cube(old_parcel_cloud->r_drop_0[p_idx]);
        //Reset radius 
        old_parcel_cloud->radius[p_idx] = old_parcel_cloud->r_drop_0[p_idx];
        //Reset r_bubble and v_bubble
        old_parcel_cloud->r_bubble[p_idx] = 0.0;
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        old_parcel_cloud->recovery_time[p_idx] = current_time;
        old_parcel_cloud->recovery_count[p_idx]++;
        int recovery_count = old_parcel_cloud->recovery_count[p_idx];
        
        printf("               [RECOVERY #%d] Recovery time logged: %.6e s\n",
               recovery_count, current_time);
        printf("               [RECOVERY #%d] Recovery count: %d\n",
               recovery_count, recovery_count);
        
        return;  // Exit RPE solver
    }
    
    // Note: Recovery success check removed - parcels are converted to children and won't re-enter RPE
    
    // Diagnostic: Check if bubble is actually growing
    static int growth_check_count = 0;
    if (growth_check_count < 5 && state.R > R_before) {
        printf("[RPE_GROWTH] R: %.6e -> %.6e (+%.3e), Rdot: %.6e -> %.6e, dRdt=%.3e, dt_sub=%.3e\n",
               R_before, state.R, state.R-R_before, Rdot_before, state.Rdot, derivs.dRdt, dt_sub);
        growth_check_count++;
    }
    
    // Check if droplet has cooled below saturation temperature (but skip if in recovery mode)
    CONVERGE_precision_t T_sat_check = T_satNH3(P_amb);
    if (state.T_drop < T_sat_check && old_parcel_cloud->recovery_time[p_idx] == 0.0) {
        // Droplet subcooled - stop bubble growth
        static int subcool_count = 0;
        if (subcool_count < 3) {
            printf("[RPE_STOP] Subcooled: T_drop=%.2f K < T_sat=%.2f K, stopping bubble growth\n",
                   state.T_drop, T_sat_check);
            subcool_count++;
        }
        old_parcel_cloud->breakup_phase[p_idx] = 16;  // Subcooled T < T_sat
        old_parcel_cloud->film_flag[p_idx] = 16;
        old_parcel_cloud->r_drop_0[p_idx] = old_parcel_cloud->radius[p_idx];
        old_parcel_cloud->r_bubble[p_idx] = 0.0;
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        return;
    }
    
    // If in recovery mode and would be subcooled, log it but don't kill
    if (state.T_drop < T_sat_check && old_parcel_cloud->recovery_time[p_idx] > 0.0) {
        static int subcool_recovery_count = 0;
        if (subcool_recovery_count < 5) {
            printf("[RPE_SUBCOOL_IN_RECOVERY] p_idx=%li, T_drop=%.2f K < T_sat=%.2f K, but protected (recovery_time=%.3e s)\n",
                   p_idx, state.T_drop, T_sat_check, old_parcel_cloud->recovery_time[p_idx]);
            subcool_recovery_count++;
        }
    }
    
    // Enforce constraint: R <= Ro
    if (state.R > params.Ro) {
        static int R_limit_count = 0;
        if (R_limit_count < 3) {
            printf("[RPE_STOP] Bubble hit droplet edge: R=%.3e >= Ro=%.3e, capping growth\n",
                   state.R, params.Ro);
            R_limit_count++;
        }
        state.R = 0.95 * params.Ro;
        state.Rdot = 0.0;  // Stop growth at droplet edge
    }
    
    // Check for very small velocity (but skip if in recovery mode)
    if (state.Rdot < 1.0e-10 && old_parcel_cloud->recovery_time[p_idx] == 0.0) {
        static int small_vel_count = 0;
        if (small_vel_count < 3) {
            printf("[RPE_STOP] Bubble velocity too small: Rdot=%.3e < 1e-10, stopping\n",
                   state.Rdot);
            small_vel_count++;
        }
        old_parcel_cloud->breakup_phase[p_idx] = 17;  // Rdot too small
        old_parcel_cloud->film_flag[p_idx] = 17;
        old_parcel_cloud->r_drop_0[p_idx] = old_parcel_cloud->radius[p_idx];
        old_parcel_cloud->r_bubble[p_idx] = 0.0;
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        return;
    }
    
    // If in recovery mode and velocity is small, log it but don't kill
    if (state.Rdot < 1.0e-10 && old_parcel_cloud->recovery_time[p_idx] > 0.0) {
        static int small_vel_recovery_count = 0;
        if (small_vel_recovery_count < 5) {
            printf("[RPE_SMALL_VEL_IN_RECOVERY] p_idx=%li, Rdot=%.3e < 1e-10, but protected (recovery_time=%.3e s)\n",
                   p_idx, state.Rdot, old_parcel_cloud->recovery_time[p_idx]);
            small_vel_recovery_count++;
        }
    }
    
    // Update parcel cloud with new state
    old_parcel_cloud->r_bubble[p_idx] = state.R;
    old_parcel_cloud->v_bubble[p_idx] = state.Rdot;
    old_parcel_cloud->temp[p_idx] = state.T_drop;
    old_parcel_cloud->m_bubble[p_idx] = state.m_b;  // CRITICAL: Store bubble mass!
    
#if RPE_DEBUG_LOGGING
    // Lifetime-based tracking: log the oldest parcel in thermal breakup
    static FILE* rpe_debug_file = NULL;
    static CONVERGE_precision_t max_lifetime_seen = 0.0;
    static int log_call_counter = 0;
    static int first_log = 1;
    
    CONVERGE_precision_t lifetime = old_parcel_cloud->lifetime[p_idx];
    int in_thermal_breakup = (old_parcel_cloud->breakup_phase[p_idx] >= 1 && 
                              old_parcel_cloud->breakup_phase[p_idx] <= 4);
    
    // Initialize file and track oldest parcel
    if (in_thermal_breakup && lifetime > max_lifetime_seen) {
        if (rpe_debug_file == NULL) {
            rpe_debug_file = fopen(LOG_FILE_NAME, "w");
            if (rpe_debug_file != NULL) {
                fprintf(rpe_debug_file, "time,lifetime,R,Rdot,T_drop,m_b,Nu,Q_conv,mdot,Pb,P_sat,P_amb,Ro,dRdt,dRdotdt,dTdt,dmbdt,rho_l,mu_l,sigma,k_l,cp_l,L_v\n");
                printf("\n[RPE_DEBUG] Opened log file: %s\n", LOG_FILE_NAME);
            }
        }
        max_lifetime_seen = lifetime;
        if (first_log) {
            printf("[RPE_DEBUG] Tracking parcel with lifetime = %.6e s\n", lifetime);
            first_log = 0;
        }
    }
    
    // Log if this is the oldest parcel and it's time to log
    if (in_thermal_breakup && lifetime >= max_lifetime_seen && 
        log_call_counter % LOG_INTERVAL == 0 && rpe_debug_file != NULL) {
        
        CONVERGE_precision_t time = CONVERGE_simulation_time();
        
        // Calculate Pb for logging
        CONVERGE_precision_t Vb_log = (4.0/3.0) * PI * state.R * state.R * state.R;
        CONVERGE_precision_t rho_v_log = safe_divide(state.m_b, Vb_log, 1e-6);
        if (rho_v_log < 1e-6) rho_v_log = 1e-6;
        CONVERGE_precision_t Pb_log = rho_v_log * params.R_spec * state.T_drop;
        
        fprintf(rpe_debug_file, 
            "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
            time, lifetime, state.R, state.Rdot, state.T_drop, state.m_b,
            Nu, Q_conv, mdot, Pb_log, P_sat, params.P_amb, params.Ro,
            derivs.dRdt, derivs.dRdotdt, derivs.dTdt, derivs.dmbdt,
            params.rho_l, params.mu_l, params.sigma, params.k_l, params.cp_l, params.L_v);
        fflush(rpe_debug_file);
    }
    log_call_counter++;
#endif
}
