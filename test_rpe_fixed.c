// Simplified test - just core physics
#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846

int main() {
    // NH3 at 355 K, 2 bar - proper superheat
    double T = 355.0;
    double P_amb = 2.0e5;
    double P_sat = 2.42e5;  // At 355 K
    double Ro = 82.5e-6;
    
    double R = 1e-6;  // Start at 1 micron (more realistic)
    double Rdot = 0.0;
    double dt = 1e-10;
    
    FILE* fp = fopen("simple_test.csv", "w");
    fprintf(fp, "time_us,R_um,Rdot_m_s\n");
    
    printf("Simple test: T=%.1f K, P_sat=%.2f bar, P_amb=%.2f bar\n", T, P_sat*1e-5, P_amb*1e-5);
    printf("Superheat = %.2f bar\n\n", (P_sat-P_amb)*1e-5);
    
    for (int i = 0; i < 50000; i++) {
        double rho_l = 610.0;
        double sigma = 0.020;
        double mu_l = 1.2e-4;
        
        // Simplified RP equation (inertial limit)
        double driving = (P_sat - P_amb) / (rho_l * R);
        double Rddot = driving - 1.5*Rdot*Rdot/R;
        
        Rdot += Rddot * dt;
        R += Rdot * dt;
        
        if (i % 500 == 0) {
            double t_us = i * dt * 1e6;
            fprintf(fp, "%.3f,%.3f,%.3f\n", t_us, R*1e6, Rdot);
            if (i % 5000 == 0) {
                printf("t=%.2f µs: R=%.3f µm, Rdot=%.2f m/s\n", t_us, R*1e6, Rdot);
            }
        }
        
        if (R > 0.9*Ro) {
            printf("\nBubble filled droplet at t=%.2f µs\n", i*dt*1e6);
            break;
        }
    }
    
    fclose(fp);
    printf("\nOutput: simple_test.csv\n");
    return 0;
}
