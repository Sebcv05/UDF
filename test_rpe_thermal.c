/*
 * RPE Unit Test - THERMAL LIMITING with Film Thickness
 * Matches Python euler_explicit.py implementation
 * 
 * Compile: gcc -o test_rpe_thermal test_rpe_thermal.c -lm
 * Run: ./test_rpe_thermal
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define R_SPEC_NH3 488.2  // J/(kg·K)

// Antoine coefficients (Stull 1947)
double P_sat_NH3(double T) {
    double A1 = 4.86886, B1 = 1113.928, C1 = -10.409;
    double A2 = 3.18757, B2 = 506.713, C2 = -80.78;
    double P_bar = pow(10.0, A1 - B1/(T + C1));
    if (P_bar < 0.9933080) P_bar = pow(10.0, A2 - B2/(T + C2));
    return P_bar * 1e5;
}

double T_sat_NH3(double P) {
    double A1 = 4.86886, B1 = 1113.928, C1 = -10.409;
    double A2 = 3.18757, B2 = 506.713, C2 = -80.78;
    double P_bar = P / 1e5;
    
    if (P_bar > 0.9933080) {
        return B1 / (A1 - log10(P_bar)) - C1;
    } else {
        return B2 / (A2 - log10(P_bar)) - C2;
    }
}

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

// Compute heat and mass transfer - DUAL GEOMETRY (matches Python)
void compute_thermal_mass_transfer(
    double R, double Rdot, double T_drop, double m_b,
    double Ro, double rho_l, double mu_l, double k_l, double cp_l, 
    double R_spec, double L_v, double max_Nu,
    double* Nu_out, double* Q_out, double* mdot_out, int* geom_out)
{
    // Calculate bubble pressure from current mass
    double Vb = (4.0/3.0) * PI * R*R*R;
    double rho_v = (Vb > 1e-30) ? (m_b / Vb) : 1e-6;
    if (rho_v < 1e-6) rho_v = 1e-6;
    
    double T_drop_safe = (T_drop > 1e-3) ? T_drop : 1e-3;
    double Pb = rho_v * R_spec * T_drop_safe;
    
    // Interface temperature from bubble pressure
    double T_int = T_sat_NH3(Pb);
    
    // Temperature difference
    double dT = T_drop - T_int;
    
    // Reynolds and Prandtl numbers
    double L_char = 2.0 * R;
    double Re = rho_l * fabs(Rdot) * L_char / mu_l;
    double Pr = mu_l * cp_l / k_l;
    
    // Nusselt number (Ranz-Marshall)
    double Nu_uncapped = 2.0 + 0.6 * sqrt(Re) * pow(Pr, 1.0/3.0);
    double Nu = (Nu_uncapped < max_Nu) ? Nu_uncapped : max_Nu;
    
    double A_bubble = 4.0 * PI * R * R;
    
    // ========== DUAL GEOMETRY ==========
    // Semi-infinite formulation
    double h_conv_si = Nu * k_l / L_char;
    double Q_si = h_conv_si * A_bubble * dT;
    
    // Film formulation
    double film_thickness = Ro - R;
    double h_conv_film, Q_film;
    
    if (film_thickness > 1e-12) {
        h_conv_film = Nu * k_l / film_thickness;
        Q_film = h_conv_film * A_bubble * dT;
    } else {
        h_conv_film = 1e20;
        Q_film = 0.0;
    }
    
    // Use maximum (least restrictive)
    double Q_conv = (Q_film > Q_si) ? Q_film : Q_si;
    int geometry = (Q_film > Q_si) ? 1 : 0;  // 0=semi-infinite, 1=film
    
    // Mass transfer rate (thermal limiting)
    double mdot = (L_v > 0) ? (Q_conv / L_v) : 0.0;
    if (mdot < 0.0) mdot = 0.0;  // Evaporation only
    
    *Nu_out = Nu;
    *Q_out = Q_conv;
    *mdot_out = mdot;
    *geom_out = geometry;
}

// Compute derivatives
void compute_derivatives(
    const BubbleState* state,
    double Ro, double rho_l, double mu_l, double k_l, double cp_l,
    double sigma, double R_spec, double L_v, double P_amb, double m_drop,
    double max_Nu, double mdot,
    BubbleDerivatives* derivs)
{
    double R = state->R;
    double Rdot = state->Rdot;
    double T_drop = state->T_drop;
    double m_b = state->m_b;
    
    if (R < 1e-12) R = 1e-12;
    
    // Calculate bubble pressure
    double Vb = (4.0/3.0) * PI * R*R*R;
    double rho_v = (Vb > 1e-30) ? (m_b / Vb) : 1e-6;
    if (rho_v < 1e-6) rho_v = 1e-6;
    
    double T_drop_safe = (T_drop > 1e-3) ? T_drop : 1e-3;
    double Pb = rho_v * R_spec * T_drop_safe;
    
    derivs->dRdt = Rdot;
    
    // Rayleigh-Plesset equation with inertia
    double pressure_term = (Pb - P_amb - 2.0*sigma/R - 4.0*mu_l*Rdot/R) / (rho_l * R);
    double inertial_term = -1.5 * Rdot * Rdot / R;
    double Rddot = pressure_term + inertial_term;
    
    if (Rddot > 1e10) Rddot = 1e10;
    if (Rddot < -1e10) Rddot = -1e10;
    
    derivs->dRdotdt = Rddot;
    derivs->dmbdt = mdot;
    
    // Temperature evolution (energy balance)
    double Nu, Q_conv, mdot_local;
    int geom;
    compute_thermal_mass_transfer(R, Rdot, T_drop, m_b, Ro, rho_l, mu_l, k_l, cp_l,
                                   R_spec, L_v, max_Nu, &Nu, &Q_conv, &mdot_local, &geom);
    
    double Q_evap = L_v * mdot;
    double dTdt = (Q_conv - Q_evap) / (m_drop * cp_l);
    
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
    printf("  RPE Unit Test - THERMAL LIMITING with Film Thickness\n");
    printf("  Matches Python euler_explicit.py\n");
    printf("=================================================================\n\n");
    
    // Parameters
    double T0 = 293.15;          // K
    double P_amb = 2.0e5;        // Pa
    double Ro = 82.5e-6;         // m
    double rho_l = 610.0;        // kg/m³
    double mu_l = 1.2e-4;        // Pa·s
    double k_l = 0.5;            // W/(m·K)
    double cp_l = 4800.0;        // J/(kg·K)
    double sigma = 0.020;        // N/m
    double L_v = 1.2e6;          // J/kg
    double max_Nu = 1000.0;
    
    double m_drop = (4.0/3.0) * PI * Ro*Ro*Ro * rho_l;
    
    // Initial state
    BubbleState state;
    state.R = 1e-6;              // Start at 1 µm
    state.Rdot = 0.0;
    state.T_drop = T0;
    
    double P_sat = P_sat_NH3(T0);
    double rho_b = P_sat / (R_SPEC_NH3 * T0);
    double Vb = (4.0/3.0) * PI * state.R*state.R*state.R;
    state.m_b = rho_b * Vb;
    
    printf("Initial Conditions:\n");
    printf("  T = %.2f K, P_sat = %.2f bar, P_amb = %.2f bar\n", T0, P_sat*1e-5, P_amb*1e-5);
    printf("  Superheat = %.2f bar\n", (P_sat - P_amb)*1e-5);
    printf("  Ro = %.2f µm, R0 = %.2f µm\n\n", Ro*1e6, state.R*1e6);
    
    FILE* fp = fopen("thermal_test.csv", "w");
    fprintf(fp, "time_us,R_um,Rdot_m_s,T_drop_K,film_um,geometry,Nu,Q_conv_W\n");
    
    double dt = 1e-10;
    double t = 0.0;
    int n_steps = 200000;  // 20 µs
    int print_interval = 1000;
    
    printf("%-8s %-10s %-10s %-8s %-10s %-6s\n", 
           "t(µs)", "R(µm)", "Rdot(m/s)", "T(K)", "film(µm)", "Geom");
    printf("-------------------------------------------------------------\n");
    
    for (int i = 0; i < n_steps; i++) {
        // Compute mass transfer
        double Nu, Q_conv, mdot;
        int geometry;
        compute_thermal_mass_transfer(state.R, state.Rdot, state.T_drop, state.m_b,
                                       Ro, rho_l, mu_l, k_l, cp_l, R_SPEC_NH3, 
                                       L_v, max_Nu, &Nu, &Q_conv, &mdot, &geometry);
        
        // Compute derivatives
        BubbleDerivatives derivs;
        compute_derivatives(&state, Ro, rho_l, mu_l, k_l, cp_l, sigma, 
                           R_SPEC_NH3, L_v, P_amb, m_drop, max_Nu, mdot, &derivs);
        
        // Write output
        if (i % print_interval == 0) {
            double film = (Ro - state.R) * 1e6;
            fprintf(fp, "%.3f,%.3f,%.3f,%.3f,%.3f,%d,%.2f,%.6e\n",
                    t*1e6, state.R*1e6, state.Rdot, state.T_drop, film, 
                    geometry, Nu, Q_conv);
            
            if (i % 10000 == 0) {
                printf("%-8.2f %-10.3f %-10.3f %-8.2f %-10.3f %-6s\n",
                       t*1e6, state.R*1e6, state.Rdot, state.T_drop, film,
                       geometry ? "Film" : "Semi");
            }
        }
        
        // Euler step
        euler_step(&state, &derivs, dt);
        t += dt;
        
        // Stop conditions
        if (state.R > 0.95 * Ro) {
            printf("\n>>> Bubble filled droplet at t = %.2f µs <<<\n", t*1e6);
            fprintf(fp, "%.3f,%.3f,%.3f,%.3f,%.3f,%d,%.2f,%.6e\n",
                    t*1e6, state.R*1e6, state.Rdot, state.T_drop, 
                    (Ro-state.R)*1e6, geometry, Nu, Q_conv);
            break;
        }
        
        double T_sat_check = T_sat_NH3(P_amb);
        if (state.T_drop < T_sat_check) {
            printf("\n>>> Subcooled at t = %.2f µs <<<\n", t*1e6);
            break;
        }
    }
    
    fclose(fp);
    
    printf("\n=================================================================\n");
    printf("Final State:\n");
    printf("  R = %.3f µm, Rdot = %.3f m/s\n", state.R*1e6, state.Rdot);
    printf("  T_drop = %.2f K (ΔT = %.2f K)\n", state.T_drop, T0 - state.T_drop);
    printf("  Time = %.2f µs\n", t*1e6);
    printf("\nOutput: thermal_test.csv\n");
    printf("=================================================================\n");
    
    return 0;
}
