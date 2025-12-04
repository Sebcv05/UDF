/*
 * test_song_rpe.c
 * 
 * Standalone test of Song et al. RPE model for isothermal bubble growth
 * Mimics song_temp_sweep.py functionality
 * 
 * Generates CSV output for plotting temperature sweep at P=2bar
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Physical constants
#define PI 3.14159265358979323846

// Ammonia liquid properties
#define RHO_L 609.6              // kg/m³
#define MU_L 138e-6              // Pa·s
#define SIGMA 0.021              // N/m
#define CP_L 4700.0              // J/(kg·K)
#define L_V 1.37e6               // J/kg
#define K_LIQ 0.54               // W/(m·K)

// Molecular properties
#define M_MOL 0.017031           // kg/mol
#define R_UNIV 8.314462618       // J/(mol·K)
#define R_SPEC (R_UNIV / M_MOL)  // 488.2 J/(kg·K)

// Antoine coefficients for ammonia (T in K, P in bar)
#define A_ANT 4.86886
#define B_ANT 1113.928
#define C_ANT (-10.409)

// Model parameters
#define KAPPA 0.0                // Surface viscosity (N/m)
#define P_R0 1.0e6               // Initial residual gas pressure (Pa)
#define MAX_NU 100.0

// Simulation parameters
#define R_DROP_0 50.0e-6         // Initial droplet radius (m)
#define TMAX 100e-6              // Max simulation time (s)
#define VOID_TARGET 0.99         // Target void fraction
#define R_BUBBLE_MULTIPLIER 1.01 // R_bubble_0 = 1.01 * R_c

// Numerical parameters
#define DT_INITIAL 1e-9          // Initial timestep (s) - increased from 1e-10
#define DT_MIN 1e-10
#define DT_MAX 1e-7
#define SAFETY_FACTOR 0.1        // More conservative - reduced from 0.9
#define TOLERANCE 1e-6

// Temperature sweep conditions (at P=2bar)
typedef struct {
    double T0;          // Temperature (K)
    double P_ambient;   // Ambient pressure (Pa)
    const char* label;  // Label for output
} Condition;

// Bubble state
typedef struct {
    double R;           // Bubble radius (m)
    double Rdot;        // Wall velocity (m/s)
    double T_drop;      // Droplet temperature (K) - isothermal
    double R0;          // Initial bubble radius (m)
    double Ro;          // Initial droplet radius (m)
} BubbleState;

// Simulation parameters
typedef struct {
    double rho_l;       // Liquid density
    double mu_l;        // Liquid viscosity
    double sigma;       // Surface tension
    double R_spec;      // Specific gas constant
    double P_amb;       // Ambient pressure
    double P_r0;        // Initial residual pressure
    double kappa;       // Surface viscosity
} SongParams;

// Result storage
typedef struct {
    double* t;          // Time array
    double* R;          // Radius array
    double* Rdot;       // Velocity array
    double* epsilon;    // Void fraction array
    double* rho_m;      // Mixture density array
    double* P_sat;      // Saturation pressure array
    int n_points;       // Number of data points
    int capacity;       // Array capacity
} ResultArrays;

// ============================================================================
// Antoine Equation Functions
// ============================================================================

double Psat_from_Antoine(double T) {
    // For ammonia: log10(P_bar) = A - B/(T_K + C)
    double log10P_bar = A_ANT - B_ANT / (T + C_ANT);
    double P_bar = pow(10.0, log10P_bar);
    return P_bar * 1e5;  // Convert bar to Pa
}

double Tsat_from_Antoine(double P) {
    // Inverse: T_K = B/(A - log10(P_bar)) - C
    double P_bar = P / 1e5;
    if (P_bar < 1e-9) P_bar = 1e-9;
    double log10P_bar = log10(P_bar);
    double denom = A_ANT - log10P_bar;
    if (fabs(denom) < 1e-9) denom = 1e-9;
    return B_ANT / denom - C_ANT;
}

// ============================================================================
// Critical Radius Calculation
// ============================================================================

double compute_critical_radius(double T, double P_amb) {
    double P_sat = Psat_from_Antoine(T);
    double dP = P_sat - P_amb;
    if (dP <= 0.0) {
        printf("WARNING: P_sat <= P_amb, cannot compute R_c\n");
        return 1e-9;
    }
    return 2.0 * SIGMA / dP;
}

// ============================================================================
// Void Fraction and Mixture Density
// ============================================================================

double compute_void_fraction(double R, double Ro) {
    double epsilon = (R*R*R) / (Ro*Ro*Ro);
    return (epsilon > 0.999) ? 0.999 : epsilon;
}

double compute_mixture_density(double epsilon, double rho_v, double rho_l) {
    return epsilon * rho_v + (1.0 - epsilon) * rho_l;
}

// ============================================================================
// Song RPE Acceleration
// ============================================================================

double compute_song_acceleration(
    const BubbleState* state,
    const SongParams* params,
    double rho_m,
    double P_sat
) {
    double R = state->R;
    double Rdot = state->Rdot;
    double R0 = state->R0;
    
    // Prevent division by zero
    if (R < 1e-12) R = 1e-12;
    
    // Initial pressure term: (2σ/R0 + P_r0) * (R0/R)^3
    double P_init = (2.0*params->sigma/R0 + params->P_r0) * pow(R0/R, 3.0);
    
    // Pressure terms: P_sat - P_∞ + P_init - 2σ/R
    double pressure_term = P_sat - params->P_amb + P_init - 2.0*params->sigma/R;
    
    // Viscous dissipation: -4μ·Ṙ/R - 4κ·Ṙ/R²
    double viscous_term = -4.0*params->mu_l*Rdot/R - 4.0*params->kappa*Rdot/(R*R);
    
    // Inertial term: -(3/2)·Ṙ²/R
    double inertial_term = -1.5 * Rdot * Rdot / R;
    
    // R̈ = [pressure + viscous] / (ρ_m·R) + inertial
    double Rddot = (pressure_term + viscous_term) / (rho_m * R) + inertial_term;
    
    return Rddot;
}

// ============================================================================
// Result Array Management
// ============================================================================

ResultArrays* create_result_arrays(int initial_capacity) {
    ResultArrays* res = (ResultArrays*)malloc(sizeof(ResultArrays));
    res->capacity = initial_capacity;
    res->n_points = 0;
    res->t = (double*)malloc(initial_capacity * sizeof(double));
    res->R = (double*)malloc(initial_capacity * sizeof(double));
    res->Rdot = (double*)malloc(initial_capacity * sizeof(double));
    res->epsilon = (double*)malloc(initial_capacity * sizeof(double));
    res->rho_m = (double*)malloc(initial_capacity * sizeof(double));
    res->P_sat = (double*)malloc(initial_capacity * sizeof(double));
    return res;
}

void resize_result_arrays(ResultArrays* res) {
    res->capacity *= 2;
    res->t = (double*)realloc(res->t, res->capacity * sizeof(double));
    res->R = (double*)realloc(res->R, res->capacity * sizeof(double));
    res->Rdot = (double*)realloc(res->Rdot, res->capacity * sizeof(double));
    res->epsilon = (double*)realloc(res->epsilon, res->capacity * sizeof(double));
    res->rho_m = (double*)realloc(res->rho_m, res->capacity * sizeof(double));
    res->P_sat = (double*)realloc(res->P_sat, res->capacity * sizeof(double));
}

void append_result(ResultArrays* res, double t, double R, double Rdot, 
                   double epsilon, double rho_m, double P_sat) {
    if (res->n_points >= res->capacity) {
        resize_result_arrays(res);
    }
    int i = res->n_points;
    res->t[i] = t;
    res->R[i] = R;
    res->Rdot[i] = Rdot;
    res->epsilon[i] = epsilon;
    res->rho_m[i] = rho_m;
    res->P_sat[i] = P_sat;
    res->n_points++;
}

void free_result_arrays(ResultArrays* res) {
    free(res->t);
    free(res->R);
    free(res->Rdot);
    free(res->epsilon);
    free(res->rho_m);
    free(res->P_sat);
    free(res);
}

// ============================================================================
// Song RPE Solver
// ============================================================================

ResultArrays* solve_song_rpe(const Condition* cond) {
    printf("\n========================================\n");
    printf("Solving: %s (T=%.2f K, P=%.2f bar)\n", 
           cond->label, cond->T0, cond->P_ambient/1e5);
    printf("========================================\n");
    
    // Initialize parameters
    SongParams params;
    params.rho_l = RHO_L;
    params.mu_l = MU_L;
    params.sigma = SIGMA;
    params.R_spec = R_SPEC;
    params.P_amb = cond->P_ambient;
    params.P_r0 = P_R0;
    params.kappa = KAPPA;
    
    // Calculate saturation pressure (isothermal)
    double P_sat = Psat_from_Antoine(cond->T0);
    double T_sat = Tsat_from_Antoine(cond->P_ambient);
    
    printf("P_sat = %.3e Pa (%.2f bar)\n", P_sat, P_sat/1e5);
    printf("T_sat = %.2f K\n", T_sat);
    printf("Superheat = %.2f K\n", cond->T0 - T_sat);
    
    // Check superheat condition
    if (P_sat <= params.P_amb) {
        printf("ERROR: Not superheated! P_sat <= P_amb\n");
        return NULL;
    }
    
    // Calculate critical radius and initial bubble radius
    double R_c = compute_critical_radius(cond->T0, params.P_amb);
    double R_bubble_0 = R_BUBBLE_MULTIPLIER * R_c;
    
    printf("R_c = %.3e m (%.3f μm)\n", R_c, R_c*1e6);
    printf("R_bubble_0 = %.3e m (%.3f μm)\n", R_bubble_0, R_bubble_0*1e6);
    
    // Initialize bubble state
    BubbleState state;
    state.R = R_bubble_0;
    state.Rdot = 0.0;
    state.T_drop = cond->T0;  // Isothermal
    state.R0 = R_bubble_0;
    state.Ro = R_DROP_0;
    
    // Calculate initial vapor density
    double rho_v = P_sat / (params.R_spec * cond->T0);
    printf("rho_v = %.3f kg/m³\n", rho_v);
    
    // Initialize result storage
    ResultArrays* results = create_result_arrays(10000);
    
    // Time integration variables
    double t = 0.0;
    double dt = DT_INITIAL;
    int step = 0;
    int log_interval = 1000;
    
    // Store initial state
    double epsilon = compute_void_fraction(state.R, state.Ro);
    double rho_m = compute_mixture_density(epsilon, rho_v, params.rho_l);
    append_result(results, t, state.R, state.Rdot, epsilon, rho_m, P_sat);
    
    printf("\nStarting integration...\n");
    printf("Initial: ε=%.6f, ρ_m=%.2f kg/m³\n", epsilon, rho_m);
    
    // Main integration loop
    while (t < TMAX) {
        // Compute current void fraction and mixture density
        epsilon = compute_void_fraction(state.R, state.Ro);
        rho_m = compute_mixture_density(epsilon, rho_v, params.rho_l);
        
        // Check termination criterion FIRST
        if (epsilon >= VOID_TARGET) {
            printf("\nReached void fraction target: ε=%.4f at t=%.3e s\n", 
                   epsilon, t);
            break;
        }
        
        // Check for collapse
        if (state.Rdot < 0.0) {
            printf("\nBubble collapse detected at t=%.3e s, Rdot=%.3e m/s\n", 
                   t, state.Rdot);
            break;
        }
        
        // Compute acceleration
        double Rddot = compute_song_acceleration(&state, &params, rho_m, P_sat);
        
        // Explicit Euler integration
        state.Rdot += Rddot * dt;
        state.R += state.Rdot * dt;
        t += dt;
        step++;
        
        // Safety checks
        if (state.R < 1e-12) state.R = 1e-12;
        
        // Cap radius at droplet boundary (but shouldn't happen if void fraction check works)
        if (state.R > 0.999 * state.Ro) {
            state.R = 0.999 * state.Ro;
            state.Rdot = 0.0;
            printf("\nBubble capped at droplet edge at t=%.3e s\n", t);
            break;
        }
        
        // Store result
        append_result(results, t, state.R, state.Rdot, epsilon, rho_m, P_sat);
        
        // Adaptive timestep (CFL-like condition based on velocity and acceleration)
        if (fabs(state.Rdot) > 1e-10) {
            double dt_CFL = SAFETY_FACTOR * state.R / fabs(state.Rdot);
            dt = fmin(fmax(dt_CFL, DT_MIN), DT_MAX);
        } else {
            // If velocity is very small, use acceleration-based timestep
            if (fabs(Rddot) > 1e-10) {
                double dt_accel = sqrt(SAFETY_FACTOR * state.R / fabs(Rddot));
                dt = fmin(fmax(dt_accel, DT_MIN), DT_MAX);
            }
        }
        
        // Periodic logging
        if (step % log_interval == 0) {
            printf("Step %d: t=%.3e s, R=%.3f μm, Rdot=%.3f m/s, ε=%.4f, dt=%.3e s\n",
                   step, t, state.R*1e6, state.Rdot, epsilon, dt);
        }
    }
    
    printf("\nSimulation complete:\n");
    printf("  Final time: %.3e s (%.2f μs)\n", t, t*1e6);
    printf("  Final R: %.3e m (%.2f μm)\n", state.R, state.R*1e6);
    printf("  Final Rdot: %.3f m/s\n", state.Rdot);
    printf("  Final ε: %.4f\n", epsilon);
    printf("  Total steps: %d\n", step);
    printf("  Data points: %d\n", results->n_points);
    
    return results;
}

// ============================================================================
// CSV Output
// ============================================================================

void write_csv_output(const char* filename, Condition* conditions, 
                      ResultArrays** all_results, int n_conditions) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        printf("ERROR: Could not open %s for writing\n", filename);
        return;
    }
    
    // Write header
    fprintf(fp, "condition,time_s,time_us,R_m,R_um,Rdot_ms,epsilon,rho_m,P_sat_Pa,P_sat_kPa,film_um\n");
    
    // Write data for each condition
    for (int c = 0; c < n_conditions; c++) {
        ResultArrays* res = all_results[c];
        if (res == NULL) continue;
        
        for (int i = 0; i < res->n_points; i++) {
            double film_um = (R_DROP_0 - res->R[i]) * 1e6;
            fprintf(fp, "%s,%.12e,%.6e,%.12e,%.6f,%.6f,%.9f,%.6f,%.6e,%.3f,%.6f\n",
                    conditions[c].label,
                    res->t[i],
                    res->t[i] * 1e6,
                    res->R[i],
                    res->R[i] * 1e6,
                    res->Rdot[i],
                    res->epsilon[i],
                    res->rho_m[i],
                    res->P_sat[i],
                    res->P_sat[i] / 1e3,
                    film_um);
        }
    }
    
    fclose(fp);
    printf("\nCSV output written to %s\n", filename);
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    printf("========================================\n");
    printf("Song et al. RPE Model - Temperature Sweep\n");
    printf("========================================\n");
    printf("\nParameters:\n");
    printf("  R_DROP_0 = %.1f μm\n", R_DROP_0*1e6);
    printf("  R_BUBBLE_0 = %.2f × R_c\n", R_BUBBLE_MULTIPLIER);
    printf("  Ambient pressure = 2.0 bar\n");
    printf("  Target void fraction = %.2f\n", VOID_TARGET);
    printf("  Max simulation time = %.1f μs\n", TMAX*1e6);
    
    // Validate Antoine coefficients
    double P_test = Psat_from_Antoine(293.15);
    printf("\nAntoine validation: P_sat(293.15K) = %.2f bar\n", P_test/1e5);
    if (P_test/1e5 < 8.0 || P_test/1e5 > 8.6) {
        printf("WARNING: Antoine coefficients may be incorrect!\n");
    }
    
    // Define conditions
    Condition conditions[] = {
        {255.15, 2e5, "T=255K"},
        {263.15, 2e5, "T=263K"},
        {273.15, 2e5, "T=273K"},
        {283.15, 2e5, "T=283K"},
        {293.15, 2e5, "T=293K"},
        {303.15, 2e5, "T=303K"},
        {313.15, 2e5, "T=313K"},
        {323.15, 2e5, "T=323K"}
    };
    int n_conditions = sizeof(conditions) / sizeof(Condition);
    
    // Allocate results storage
    ResultArrays** all_results = (ResultArrays**)malloc(n_conditions * sizeof(ResultArrays*));
    
    // Solve for each condition
    for (int i = 0; i < n_conditions; i++) {
        all_results[i] = solve_song_rpe(&conditions[i]);
    }
    
    // Write CSV output
    write_csv_output("song_temp_sweep_c.csv", conditions, all_results, n_conditions);
    
    // Free memory
    for (int i = 0; i < n_conditions; i++) {
        if (all_results[i]) free_result_arrays(all_results[i]);
    }
    free(all_results);
    
    printf("\nDone!\n");
    printf("Use Python script to plot: plot_song_c_results.py\n");
    
    return 0;
}
