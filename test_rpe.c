/*
 * RPE_euler Unit Test
 * 
 * Standalone test of RPE solver with fixed inputs
 * Tests bubble growth from nucleation to near-droplet size
 * 
 * Compile: gcc -o test_rpe test_rpe.c -lm
 * Run: ./test_rpe
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
    // NIST coefficients for ammonia (K, Pa)
    double A = 9.96268;
    double B = 1617.9;
    double C = 6.65;
    return exp(A - B / (T + C)) * 1000.0;
}

double T_sat_NH3(double P) {
    // Inverse of Antoine equation
    double A = 9.96268;
    double B = 1617.9;
    double C = 6.65;
    return B / (A - log(P / 1000.0)) - C;
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
    
    double h_conv = Nu * params->k_l / L_char;
    double A_bubble = 4.0 * PI * R * R;
    double Q_conv = h_conv * A_bubble * dT;
    
    double mdot = Q_conv / params->L_v;
    
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
    printf("  RPE_euler Unit Test - NH3 Bubble Growth\n");
    printf("=================================================================\n\n");
    
    // Test case: NH3 droplet at 323 K, 2 bar ambient
    RPE_Params params;
    params.rho_l = 610.0;        // kg/m³ (NH3 liquid at 323K)
    params.mu_l = 1.2e-4;        // Pa·s
    params.k_l = 0.5;            // W/(m·K)
    params.cp_l = 4800.0;        // J/(kg·K)
    params.sigma = 0.020;        // N/m
    params.R_spec = R_SPEC_NH3;
    params.L_v = 1.2e6;          // J/kg
    params.Ro = 82.5e-6;         // m (82.5 µm)
    params.P_amb = 2.0e5;        // Pa (2 bar)
    params.T_amb = 298.0;        // K
    params.max_Nu = 1000.0;
    params.m_drop = (4.0/3.0) * PI * params.Ro * params.Ro * params.Ro * params.rho_l;
    
    // Initial state
    BubbleState state;
    state.R = 1.0e-9;            // m (1 nm - nucleation size)
    state.Rdot = 0.0;            // m/s
    state.T_drop = 355.0;        // K (superheated)
    
    double P_sat = P_sat_NH3(state.T_drop);
    double rho_b = bubble_density_NH3(P_sat, state.T_drop);
    double Vb = (4.0/3.0) * PI * state.R * state.R * state.R;
    state.m_b = rho_b * Vb;
    
    printf("Initial Conditions:\n");
    printf("  Droplet: Ro = %.3f µm, T = %.2f K, m = %.3e kg\n", 
           params.Ro*1e6, state.T_drop, params.m_drop);
    printf("  Bubble:  R = %.3f nm, Rdot = %.3f m/s\n", state.R*1e9, state.Rdot);
    printf("  Ambient: P = %.2f bar, T = %.2f K\n", params.P_amb*1e-5, params.T_amb);
    printf("  P_sat = %.2f bar (superheat = %.2f bar)\n\n", 
           P_sat*1e-5, (P_sat - params.P_amb)*1e-5);
    
    // Open output file
    FILE* fp = fopen("rpe_unit_test.csv", "w");
    fprintf(fp, "time_us,R_um,Rdot_m_s,T_drop_K,m_b_ng,Nu,Q_conv_W,mdot_mg_s,Pb_bar,P_sat_bar,dRdt,dRdotdt,dTdt\n");
    
    // Time integration
    double dt = 1e-10;           // 0.1 ns timestep
    double t = 0.0;
    int n_steps = 100000;        // 10 µs total
    int print_interval = 1000;   // Print every 1 µs
    
    printf("Time Integration (dt = %.2e s, t_end = %.2e s):\n", dt, n_steps*dt);
    printf("%-10s %-12s %-12s %-10s %-10s %-8s %-12s\n",
           "t (µs)", "R (µm)", "Rdot (m/s)", "T (K)", "Nu", "kb*", "Status");
    printf("------------------------------------------------------------------------------------\n");
    
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
        
        // Simplified kb estimate (eta / liquid shell thickness)
        double eta_0 = 0.05 * params.Ro;  // Typical initial perturbation
        double kb_approx = (state.R > 0.01*params.Ro) ? 
                           (eta_0 / (params.Ro - state.R)) : 0.0;
        
        // Write to file
        fprintf(fp, "%.6f,%.6f,%.6f,%.6f,%.6f,%.3f,%.6e,%.6e,%.6f,%.6f,%.6e,%.6e,%.6e\n",
                t*1e6, state.R*1e6, state.Rdot, state.T_drop, state.m_b*1e12,
                Nu, Q_conv, mdot*1e6, Pb*1e-5, P_sat_current*1e-5,
                derivs.dRdt, derivs.dRdotdt, derivs.dTdt);
        
        // Print progress
        if (i % print_interval == 0) {
            char status[20];
            if (state.R > 0.9 * params.Ro) {
                sprintf(status, "FILLED");
            } else if (kb_approx > 1.0) {
                sprintf(status, "BREAKUP");
            } else {
                sprintf(status, "GROWING");
            }
            
            printf("%-10.2f %-12.3f %-12.3f %-10.2f %-10.2f %-8.3f %-12s\n",
                   t*1e6, state.R*1e6, state.Rdot, state.T_drop, Nu, kb_approx, status);
        }
        
        // Euler step
        euler_step(&state, &derivs, dt);
        t += dt;
        
        // Stop conditions
        if (state.R > 0.95 * params.Ro) {
            printf("\n>>> Bubble filled droplet at t = %.3f µs <<<\n", t*1e6);
            break;
        }
        
        double T_sat_check = T_sat_NH3(params.P_amb);
        if (state.T_drop < T_sat_check) {
            printf("\n>>> Droplet subcooled at t = %.3f µs (T = %.2f K < T_sat = %.2f K) <<<\n",
                   t*1e6, state.T_drop, T_sat_check);
            break;
        }
    }
    
    fclose(fp);
    
    printf("\n=================================================================\n");
    printf("Final State:\n");
    printf("  R = %.3f µm (%.1f%% of Ro)\n", state.R*1e6, 100.0*state.R/params.Ro);
    printf("  Rdot = %.3f m/s\n", state.Rdot);
    printf("  T_drop = %.2f K (ΔT = %.2f K)\n", state.T_drop, 323.0 - state.T_drop);
    printf("  Time = %.3f µs\n", t*1e6);
    printf("\nOutput written to: rpe_unit_test.csv\n");
    printf("=================================================================\n");
    
    return 0;
}
