#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846

double P_sat_NH3(double T) {
    double A1 = 4.86886, B1 = 1113.928, C1 = -10.409;
    double A2 = 3.18757, B2 = 506.713, C2 = -80.78;
    double P_bar = pow(10.0, A1 - B1/(T + C1));
    if (P_bar < 0.9933080) P_bar = pow(10.0, A2 - B2/(T + C2));
    return P_bar * 1e5;
}

int main() {
    double T = 293.15, P_amb = 2.0e5, Ro = 82.5e-6;
    double P_sat = P_sat_NH3(T);
    double R = 1e-6, Rdot = 0.0, dt = 1e-10;
    double rho_l = 610.0, sigma = 0.020, mu_l = 1.2e-4;
    
    printf("NH3 Bubble Growth Test (CORRECT Antoine coefficients)\n");
    printf("T=%.2f K, P_sat=%.2f bar, P_amb=%.2f bar, Superheat=%.2f bar\n\n",
           T, P_sat*1e-5, P_amb*1e-5, (P_sat-P_amb)*1e-5);
    
    FILE* fp = fopen("bubble_growth.csv", "w");
    fprintf(fp, "time_us,R_um,Rdot_m_s\n");
    
    for (int i = 0; i < 100000; i++) {
        double driving = (P_sat - P_amb - 2.0*sigma/R) / (rho_l * R);
        double Rddot = driving - 1.5*Rdot*Rdot/R;
        Rdot += Rddot * dt;
        R += Rdot * dt;
        
        if (i % 1000 == 0) {
            fprintf(fp, "%.3f,%.3f,%.3f\n", i*dt*1e6, R*1e6, Rdot);
            if (i % 10000 == 0) {
                printf("t=%5.2f µs: R=%6.2f µm, Rdot=%6.2f m/s\n", i*dt*1e6, R*1e6, Rdot);
            }
        }
        
        if (R > 0.9*Ro) {
            printf("\n>>> Bubble filled droplet at t=%.2f µs <<<\n", i*dt*1e6);
            fprintf(fp, "%.3f,%.3f,%.3f\n", i*dt*1e6, R*1e6, Rdot);
            break;
        }
    }
    
    fclose(fp);
    printf("Output: bubble_growth.csv\n");
    return 0;
}
