//RPE_euler.c
// Rayleigh-Plesset Equation solver with Explicit Euler integration

#include "lagrangian/env.h"
#include <RPE_euler.h>
#include <PsatNH3.h>
#include <TsatNH3.h>
#include <BubbleDensityNH3.h>
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
    // user_lag_var_i is used as the collapse recovery counter
    // dgre_cycle_count is used as cooldown timer (only active after recovery)
    
    // Check if this parcel is in cooldown period after recovery
    if (old_parcel_cloud->dgre_cycle_count[p_idx] > 0) {
        old_parcel_cloud->dgre_cycle_count[p_idx]--;
        
        static int cooldown_print_count = 0;
        if (cooldown_print_count < 10) {
            printf("[RPE_COOLDOWN] Parcel %ld: Skipping RPE for %d more cycles (recovery cooldown)\n", 
                   p_idx, old_parcel_cloud->dgre_cycle_count[p_idx]);
            cooldown_print_count++;
        }
        
        // Don't run RPE solver during cooldown for THIS parcel only
        return;
    }
    
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
        printf("[RPE_ERROR] Droplet radius too small: Ro=%.3e m, skipping RPE solver\n", params.Ro);
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        old_parcel_cloud->pbt[p_idx] = 0;
        old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
        return;
    }
    
    // Calculate droplet mass
    params.m_drop = (4.0/3.0) * PI * params.Ro * params.Ro * params.Ro * params.rho_l;
    
    // Get temperature-dependent properties from tables
    CONVERGE_precision_t Td = old_parcel_cloud->temp[p_idx];
    
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
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        old_parcel_cloud->pbt[p_idx] = 0;
        old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
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
    
    // SAFETY: Prevent negative Rdot (bubble collapse) with recovery logic
    if (state.Rdot < 0.0) {
        // Calculate saturation pressure
        CONVERGE_precision_t P_sat_calc;
        Saturation_PressureNH3(state.T_drop, &P_sat_calc);
        CONVERGE_precision_t T_sat_calc = T_satNH3(params.P_amb);
        
        // Always print collapse diagnostics
        printf("[RPE_STOP] Negative Rdot=%.3e, attempting recovery (bubble collapsing)\n", state.Rdot);
        printf("           T_drop=%.2f K, T_sat(P_amb)=%.2f K, P_sat(T_drop)=%.3e Pa, P_amb=%.3e Pa\n",
               state.T_drop, T_sat_calc, P_sat_calc, params.P_amb);
        printf("           R=%.3e m, Ro=%.3e m, dRdt=%.3e m/s, dRdotdt=%.3e m/s²\n",
               state.R, params.Ro, derivs.dRdt, derivs.dRdotdt);
        
        // COLLAPSE RECOVERY LOGIC
        // Increment collapse counter using user_lag_var_i (clean field, not persisted weirdly)
        old_parcel_cloud->user_lag_var_i[p_idx]++;
        int collapse_count = old_parcel_cloud->user_lag_var_i[p_idx];
        
        printf("           Collapse count: %d\n", collapse_count);
        
        if (collapse_count <= 5) {
            // ATTEMPT RECOVERY: Shrink bubble to near critical and reset droplet
            CONVERGE_precision_t R_c = 2.0 * params.sigma / (P_sat_calc - params.P_amb);
            if (R_c < 1e-12) R_c = 1e-12;
            
            CONVERGE_precision_t new_R_bubble = 1.1 * R_c;
            CONVERGE_precision_t new_R_drop = old_parcel_cloud->r_drop_0[p_idx];
            
            // Calculate old and new masses
            CONVERGE_precision_t V_drop_old = (4.0/3.0) * PI * params.Ro * params.Ro * params.Ro;
            CONVERGE_precision_t V_bubble_old = (4.0/3.0) * PI * state.R * state.R * state.R;
            CONVERGE_precision_t V_liquid_old = V_drop_old - V_bubble_old;
            CONVERGE_precision_t mass_liquid_old = V_liquid_old * params.rho_l;
            
            CONVERGE_precision_t V_drop_new = (4.0/3.0) * PI * new_R_drop * new_R_drop * new_R_drop;
            CONVERGE_precision_t V_bubble_new = (4.0/3.0) * PI * new_R_bubble * new_R_bubble * new_R_bubble;
            CONVERGE_precision_t V_liquid_new = V_drop_new - V_bubble_new;
            CONVERGE_precision_t mass_per_drop_new = V_liquid_new * params.rho_l;
            
            // Conserve mass by adjusting num_drop
            CONVERGE_precision_t total_mass = old_parcel_cloud->num_drop[p_idx] * mass_liquid_old;
            CONVERGE_precision_t new_num_drop = total_mass / mass_per_drop_new;
            
            // Apply recovery
            old_parcel_cloud->r_bubble[p_idx] = new_R_bubble;
            old_parcel_cloud->radius[p_idx] = new_R_drop;
            old_parcel_cloud->num_drop[p_idx] = new_num_drop;
            old_parcel_cloud->v_bubble[p_idx] = 0.0;
            old_parcel_cloud->r_bubble_0[p_idx] = new_R_bubble;
            old_parcel_cloud->thermal_breakup_flag[p_idx] = 888;  // Signal: recovery attempted, skip rest of timestep
            
            // Set 10-cycle cooldown for THIS parcel only
            old_parcel_cloud->dgre_cycle_count[p_idx] = 10;
            
            printf("           [RECOVERY %d/5] R_bubble: %.3e -> %.3e m (1.1*Rc=%.3e)\n",
                   collapse_count, state.R, new_R_bubble, R_c);
            printf("           [RECOVERY %d/5] R_drop: %.3e -> %.3e m (r_drop_0)\n",
                   collapse_count, params.Ro, new_R_drop);
            printf("           [RECOVERY %d/5] num_drop: %.3e -> %.3e (mass conserved)\n",
                   collapse_count, old_parcel_cloud->num_drop[p_idx] * mass_liquid_old / total_mass, new_num_drop);
            printf("           [RECOVERY %d/5] Setting 10-cycle cooldown for this parcel\n",
                   collapse_count);
            
            return;  // Exit RPE solver, break out of sub-cycling loop
            
        } else {
            // GIVE UP: After 5 attempts, treat as solid droplet
            printf("           [COLLAPSE FINAL] After %d attempts, treating as solid droplet\n", collapse_count);
            
            old_parcel_cloud->r_bubble[p_idx] = 0.0;
            old_parcel_cloud->r_bubble_0[p_idx] = 0.0;
            old_parcel_cloud->radius[p_idx] = old_parcel_cloud->r_drop_0[p_idx];
            old_parcel_cloud->v_bubble[p_idx] = 0.0;
            old_parcel_cloud->pbt[p_idx] = 0;
            old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
            old_parcel_cloud->film_flag[p_idx] = 1;  // Mark as child for evaporation
            old_parcel_cloud->user_lag_var_i[p_idx] = 0;  // Reset counter
            
            return;
        }
    }
    
    // RESET collapse counter if bubble is growing successfully after recovery
    if (old_parcel_cloud->user_lag_var_i[p_idx] > 0 && state.Rdot > 0.0) {
        printf("[RPE_RECOVERY_SUCCESS] Bubble recovered! Counter %d -> 0, Rdot=%.3e m/s\n",
               old_parcel_cloud->user_lag_var_i[p_idx], state.Rdot);
        old_parcel_cloud->user_lag_var_i[p_idx] = 0;  // Reset counter on successful recovery
    }
    
    // Diagnostic: Check if bubble is actually growing
    static int growth_check_count = 0;
    if (growth_check_count < 5 && state.R > R_before) {
        printf("[RPE_GROWTH] R: %.6e -> %.6e (+%.3e), Rdot: %.6e -> %.6e, dRdt=%.3e, dt_sub=%.3e\n",
               R_before, state.R, state.R-R_before, Rdot_before, state.Rdot, derivs.dRdt, dt_sub);
        growth_check_count++;
    }
    
    // Check if droplet has cooled below saturation temperature
    CONVERGE_precision_t T_sat_check = T_satNH3(P_amb);
    if (state.T_drop < T_sat_check) {
        // Droplet subcooled - stop bubble growth
        static int subcool_count = 0;
        if (subcool_count < 3) {
            printf("[RPE_STOP] Subcooled: T_drop=%.2f K < T_sat=%.2f K, stopping bubble growth\n",
                   state.T_drop, T_sat_check);
            subcool_count++;
        }
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        old_parcel_cloud->pbt[p_idx] = 0;
        old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
        return;
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
    
    // Check for very small velocity
    if (state.Rdot < 1.0e-10) {
        static int small_vel_count = 0;
        if (small_vel_count < 3) {
            printf("[RPE_STOP] Bubble velocity too small: Rdot=%.3e < 1e-10, stopping\n",
                   state.Rdot);
            small_vel_count++;
        }
        old_parcel_cloud->v_bubble[p_idx] = 0.0;
        old_parcel_cloud->pbt[p_idx] = 0;
        old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
        return;
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
    int in_thermal_breakup = (old_parcel_cloud->pbt[p_idx] == 1);
    
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
