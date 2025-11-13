/*
 * RPE_euler Parameter Sweep Test
 * 
 * Sweeps initial droplet temperature from 273 K to 323 K (matching euler_explicit.py)
 * Fixed ambient pressure at 2 bar
 * Initial bubble radius: 1.1 * R_crit (matches CONVERGE UDF initialization)
 * Outputs Rdot vs film thickness and log(time) for each case
 * 
 * Compile: gcc -o test_rpe_sweep test_rpe.c -lm
 * Run: ./test_rpe_sweep
 * Plot: python3 plot_rpe_sweep.py
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define R_SPEC_NH3 488.2  // J/(kg·K)

// Simplified structures (matching RPE_euler.c)
typedef struct {
    double rho_l;       // Liquid density (kg/m³)
    double mu_l;        // Dynamic viscosity (Pa·s)
    double k_l;         // Thermal conductivity (W/(m·K))
    double cp_l;        // Specific heat (J/(kg·K))
    double sigma;       // Surface tension (N/m)
    double R_spec;      // Specific gas constant (J/(kg·K))
    double L_v;         // Latent heat (J/kg)
    double Ro;          // Droplet radius (m)
    double m_drop;      // Droplet mass (kg)
    double P_amb;       // Ambient pressure (Pa)
    double T_amb;       // Ambient temperature (K)
    double max_Nu;      // Maximum Nusselt number
} RPE_Params;

typedef struct {
    double R;           // Bubble radius (m)
    double Rdot;        // Bubble wall velocity (m/s)
    double T_drop;      // Droplet temperature (K)
    double m_b;         // Bubble vapor mass (kg)
} BubbleState;

typedef struct {
    double dRdt;
    double dRdotdt;
    double dTdt;
    double dmbdt;
} BubbleDerivatives;

// Simplified Antoine equation for NH3
double P_sat_NH3(double T) {
    // Antoine equation: log10(P_kPa) = A - B/(T-C)
    // Coefficients from euler_explicit.py (CORRECT for NH3)
    double A = 6.67956;
    double B = 1002.711;
    double C = 25.215;
    double log10_P_kPa = A - B / (T - C);
    return pow(10.0, log10_P_kPa) * 1000.0;  // Convert kPa to Pa
}

double T_sat_NH3(double P) {
    // Inverse of Antoine equation: T = B/(A - log10(P_kPa)) + C
    double A = 6.67956;
    double B = 1002.711;
    double C = 25.215;
    double P_kPa = P / 1000.0;  // Convert Pa to kPa
    double log10_P_kPa = log10(P_kPa);
    return B / (A - log10_P_kPa) + C;
}

// Simplified bubble density
double bubble_density_NH3(double P, double T) {
    return P / (R_SPEC_NH3 * T);
}

// Compute heat and mass transfer
void compute_thermal_mass_transfer(
    double R, double Rdot, double T_drop, double m_b,
    const RPE_Params* params,
    double* Nu_out, double* Q_out, double* mdot_out)
{
    double Vb = (4.0/3.0) * PI * R*R*R;
    double rho_v = (Vb > 1e-30) ? (m_b / Vb) : 1e-6;
    if (rho_v < 1e-6) rho_v = 1e-6;
    
    double Pb = rho_v * params->R_spec * T_drop;
    double T_int = T_sat_NH3(Pb);
    double dT = T_drop - T_int;
    
    double L_char = 2.0 * R;
    double Re = params->rho_l * fabs(Rdot) * L_char / params->mu_l;
    double Pr = params->mu_l * params->cp_l / params->k_l;
    
    double Nu = 2.0 + 0.6 * sqrt(Re) * pow(Pr, 1.0/3.0);
    if (Nu > params->max_Nu) Nu = params->max_Nu;
    
    double A_bubble = 4.0 * PI * R * R;
    
    // ========== DUAL GEOMETRY: Semi-infinite vs Film ==========
    // Match Python euler_explicit.py and CONVERGE UDF: evaluate both, use max
    
    // Semi-infinite formulation
    double h_conv_si = Nu * params->k_l / L_char;
    double Q_si = h_conv_si * A_bubble * dT;
    
    // Film formulation
    double film_thickness = params->Ro - R;
    double h_conv_film, Q_film;
    
    if (film_thickness > 1e-12) {
        h_conv_film = Nu * params->k_l / film_thickness;
        Q_film = h_conv_film * A_bubble * dT;
    } else {
        h_conv_film = 1e20;
        Q_film = 0.0;
    }
    
    // Use maximum (least restrictive) for heat transfer
    double Q_conv = (Q_film > Q_si) ? Q_film : Q_si;
    
    // Mass transfer rate (thermal limiting)
    double mdot = Q_conv / params->L_v;
    if (mdot < 0.0) mdot = 0.0;  // Evaporation only
    
    *Nu_out = Nu;
    *Q_out = Q_conv;
    *mdot_out = mdot;
}

// Compute derivatives
void compute_derivatives(
    const BubbleState* state,
    const RPE_Params* params,
    BubbleDerivatives* derivs,
    double mdot)
{
    double R = state->R;
    double Rdot = state->Rdot;
    double T_drop = state->T_drop;
    double m_b = state->m_b;
    
    if (R < 1e-12) R = 1e-12;
    
    double Vb = (4.0/3.0) * PI * R*R*R;
    double rho_v = (Vb > 1e-30) ? (m_b / Vb) : 1e-6;
    if (rho_v < 1e-6) rho_v = 1e-6;
    
    double Pb = rho_v * params->R_spec * T_drop;
    
    derivs->dRdt = Rdot;
    
    double pressure_term = (Pb - params->P_amb - 2.0*params->sigma/R - 4.0*params->mu_l*Rdot/R) 
                          / (params->rho_l * R);
    double inertial_term = -1.5 * Rdot * Rdot / R;
    double Rddot = pressure_term + inertial_term;
    
    if (Rddot > 1e10) Rddot = 1e10;
    if (Rddot < -1e10) Rddot = -1e10;
    
    derivs->dRdotdt = Rddot;
    derivs->dmbdt = mdot;
    
    double Nu, Q_conv, mdot_local;
    compute_thermal_mass_transfer(R, Rdot, T_drop, m_b, params, &Nu, &Q_conv, &mdot_local);
    
    double Q_evap = params->L_v * mdot;
    double dTdt = (Q_conv - Q_evap) / (params->m_drop * params->cp_l);
    
    if (dTdt > 1e6) dTdt = 1e6;
    if (dTdt < -1e6) dTdt = -1e6;
    
    derivs->dTdt = dTdt;
}

// Euler step
void euler_step(BubbleState* state, const BubbleDerivatives* derivs, double dt)
{
    state->R += derivs->dRdt * dt;
    state->Rdot += derivs->dRdotdt * dt;
    state->T_drop += derivs->dTdt * dt;
    state->m_b += derivs->dmbdt * dt;
    
    if (state->R < 1e-12) state->R = 1e-12;
    if (state->m_b < 0.0) state->m_b = 0.0;
    if (state->T_drop < 100.0) state->T_drop = 100.0;
    if (state->T_drop > 500.0) state->T_drop = 500.0;
}

int main() {
    printf("=================================================================\n");
    printf("  RPE_euler Parameter Sweep - NH3 Bubble Growth\n");
    printf("  P_amb = 2.0 bar (T_sat = 240 K), T0 = 273:10:323 K\n");
    printf("=================================================================\n\n");
    
    // Temperature sweep parameters (273 K to 323 K to match euler_explicit.py)
    // Note: Some cases may be subcooled at 2 bar (T_sat ~ 240 K)
    double T0_values[] = {273.0, 283.0, 293.0, 303.0, 313.0, 323.0};
    int n_temps = 6;
    double P_amb = 2.0e5;  // Pa (2 bar)
    
    // Fixed parameters
    double Ro = 82.5e-6;  // m
    
    // Open output files
    FILE* fp_all = fopen("rpe_sweep_all.csv", "w");
    fprintf(fp_all, "T0_K,time_us,R_um,Rdot_m_s,T_drop_K,film_thick_um,log10_time_us,Pb_bar,P_sat_bar\n");
    
    FILE* fp_summary = fopen("rpe_sweep_summary.csv", "w");
    fprintf(fp_summary, "T0_K,max_Rdot_m_s,time_at_max_us,R_at_max_um,film_thick_at_max_um,superheat_K\n");
    
    // Loop over temperatures
    for (int temp_idx = 0; temp_idx < n_temps; temp_idx++) {
        double T0 = T0_values[temp_idx];
        
        printf("\n-----------------------------------------------------------------\n");
        printf("Case %d/%d: T0 = %.1f K, P_amb = %.2f bar\n", 
               temp_idx+1, n_temps, T0, P_amb*1e-5);
        printf("-----------------------------------------------------------------\n");
        
        // Set up parameters (temperature-dependent properties)
        RPE_Params params;
        params.rho_l = 682.6 - 0.5 * (T0 - 273.0);  // Approximate NH3 liquid density
        params.mu_l = 1.5e-4 * exp(-0.02 * (T0 - 273.0));  // Viscosity decreases with T
        params.k_l = 0.5;
        params.cp_l = 4800.0;
        params.sigma = 0.025 - 0.0001 * (T0 - 273.0);  // Surface tension decreases with T
        params.R_spec = R_SPEC_NH3;
        params.L_v = 1.37e6 - 1000.0 * (T0 - 273.0);  // Latent heat decreases with T
        params.Ro = Ro;
        params.P_amb = P_amb;
        params.T_amb = 298.0;
        params.max_Nu = 1000.0;
        params.m_drop = (4.0/3.0) * PI * params.Ro * params.Ro * params.Ro * params.rho_l;
        
        // Calculate saturation pressure and critical radius
        double P_sat = P_sat_NH3(T0);
        
        // Critical radius calculation: Rc = 2*sigma / (P_sat - P_amb)
        double superheat_pressure = P_sat - P_amb;
        double Rc = 2.0 * params.sigma / superheat_pressure;
        
        // Initial state - start with 1.1*Rc (just above critical radius, matching UDF)
        BubbleState state;
        state.R = 1.1 * Rc;  // m (1.1 times critical radius)
        state.Rdot = 0.0;
        state.T_drop = T0;
        
        // Initialize bubble mass at equilibrium
        double rho_b = bubble_density_NH3(P_sat, state.T_drop);
        double Vb = (4.0/3.0) * PI * state.R * state.R * state.R;
        state.m_b = rho_b * Vb;
        
        double superheat = superheat_pressure;
        
        printf("  P_sat(T0) = %.3f bar, Superheat = %.3f bar\n", 
               P_sat*1e-5, superheat*1e-5);
        printf("  R_crit = %.3e m (%.3f µm), R_init = %.3e m (%.3f µm)\n",
               Rc, Rc*1e6, state.R, state.R*1e6);
        
        if (superheat < 0.0) {
            printf("  >>> SUBCOOLED: Skipping this case <<<\n");
            fprintf(fp_summary, "%.1f,0.0,0.0,0.0,0.0,%.3f\n", T0, superheat*1e-5);
            continue;
        }
        
        // Time integration - longer time to see full growth
        double dt = 1e-9;  // 1 ns (larger timestep for stability)
        double t = 0.0;
        int n_steps = 100000;  // 100 µs total
        int output_interval = 100;  // Output every 100 ns
        
        double max_Rdot = 0.0;
        double t_at_max_Rdot = 0.0;
        double R_at_max_Rdot = 0.0;
        double film_at_max_Rdot = 0.0;
        
        int step_count = 0;
        for (int i = 0; i < n_steps; i++) {
            // Compute derivatives
            BubbleDerivatives derivs;
            double Nu, Q_conv, mdot;
            compute_thermal_mass_transfer(state.R, state.Rdot, state.T_drop, state.m_b,
                                           &params, &Nu, &Q_conv, &mdot);
            compute_derivatives(&state, &params, &derivs, mdot);
            
            // Calculate diagnostic values
            double Vb_diag = (4.0/3.0) * PI * state.R * state.R * state.R;
            double rho_v_diag = (Vb_diag > 1e-30) ? (state.m_b / Vb_diag) : 1e-6;
            double Pb = rho_v_diag * params.R_spec * state.T_drop;
            double P_sat_current = P_sat_NH3(state.T_drop);
            double film_thickness = params.Ro - state.R;  // m
            
            // Track maximum Rdot
            if (fabs(state.Rdot) > max_Rdot) {
                max_Rdot = fabs(state.Rdot);
                t_at_max_Rdot = t;
                R_at_max_Rdot = state.R;
                film_at_max_Rdot = film_thickness;
            }
            
            // Output to file
            if (i % output_interval == 0) {
                double log10_time = (t > 1e-20) ? log10(t * 1e6) : -20.0;  // log10(time in µs)
                fprintf(fp_all, "%.1f,%.6e,%.6e,%.6e,%.6f,%.6e,%.6f,%.6f,%.6f\n",
                        T0, t*1e6, state.R*1e6, state.Rdot, state.T_drop,
                        film_thickness*1e6, log10_time, Pb*1e-5, P_sat_current*1e-5);
            }
            
            // Euler step
            euler_step(&state, &derivs, dt);
            t += dt;
            step_count++;
            
            // Stop conditions
            if (state.R > 0.95 * params.Ro) {
                printf("  >>> Bubble filled droplet at t = %.3f µs <<<\n", t*1e6);
                break;
            }
            
            double T_sat_check = T_sat_NH3(params.P_amb);
            if (state.T_drop < T_sat_check) {
                printf("  >>> Droplet subcooled at t = %.3f µs <<<\n", t*1e6);
                break;
            }
            
            if (t > 200e-6) {  // 200 µs timeout
                printf("  >>> Timeout at t = %.1f µs <<<\n", t*1e6);
                break;
            }
        }
        
        printf("  Final: R = %.3f µm (%.1f%%), Rdot = %.3f m/s, t = %.3f µs\n",
               state.R*1e6, 100.0*state.R/params.Ro, state.Rdot, t*1e6);
        printf("  Max Rdot: %.3f m/s at t = %.3f µs (film = %.3f µm)\n",
               max_Rdot, t_at_max_Rdot*1e6, film_at_max_Rdot*1e6);
        
        // Write summary
        fprintf(fp_summary, "%.1f,%.6e,%.6e,%.6e,%.6e,%.6f\n",
                T0, max_Rdot, t_at_max_Rdot*1e6, R_at_max_Rdot*1e6, 
                film_at_max_Rdot*1e6, superheat*1e-5);
    }
    
    fclose(fp_all);
    fclose(fp_summary);
    
    printf("\n=================================================================\n");
    printf("Parameter sweep complete!\n");
    printf("  Full data: rpe_sweep_all.csv\n");
    printf("  Summary:   rpe_sweep_summary.csv\n");
    printf("=================================================================\n");
    
    return 0;
}
