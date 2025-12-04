#!/usr/bin/env python3
"""
Plot results from test_song_rpe.c
Produces figures matching song_temp_sweep.py output
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
plt.rcParams['text.usetex'] = False  # Disable LaTeX

# Try to use science plots style (but disable LaTeX)
try:
    plt.style.use(['science', 'no-latex'])
except:
    pass

# Read CSV data
print("Reading song_temp_sweep_c.csv...")
df = pd.read_csv('song_temp_sweep_c.csv')

# Get unique conditions
conditions = df['condition'].unique()
print(f"Found {len(conditions)} conditions: {', '.join(conditions)}")

# Create main figure with 9 subplots
fig = plt.figure(figsize=(15, 12))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

ax1 = fig.add_subplot(gs[0, 0])  # R vs t
ax2 = fig.add_subplot(gs[0, 1])  # Rdot vs t
ax3 = fig.add_subplot(gs[0, 2])  # Epsilon vs t
ax4 = fig.add_subplot(gs[1, 0])  # R vs Rdot
ax5 = fig.add_subplot(gs[1, 1])  # rho_m vs t
ax6 = fig.add_subplot(gs[1, 2])  # Film thickness vs t
ax7 = fig.add_subplot(gs[2, 0])  # P_sat vs t
ax8 = fig.add_subplot(gs[2, 1])  # Epsilon vs R
ax9 = fig.add_subplot(gs[2, 2])  # Temperature (isothermal)

# Plot data for each condition
for cond in conditions:
    cond_data = df[df['condition'] == cond]
    
    # ax1: R vs time
    ax1.plot(cond_data['time_us'], cond_data['R_um'], label=cond, linewidth=1)
    
    # ax2: Rdot vs time
    ax2.plot(cond_data['time_us'], cond_data['Rdot_ms'], label=cond, linewidth=1)
    
    # ax3: Epsilon vs time
    ax3.plot(cond_data['time_us'], cond_data['epsilon'], label=cond, linewidth=1)
    
    # ax4: R vs Rdot
    ax4.plot(cond_data['R_um'], cond_data['Rdot_ms'], label=cond, linewidth=1)
    
    # ax5: rho_m vs time
    ax5.plot(cond_data['time_us'], cond_data['rho_m'], label=cond, linewidth=1)
    
    # ax6: Film thickness vs time
    ax6.plot(cond_data['time_us'], cond_data['film_um'], label=cond, linewidth=1)
    
    # ax7: P_sat vs time
    ax7.plot(cond_data['time_us'], cond_data['P_sat_kPa'], label=cond, linewidth=1)
    
    # ax8: Epsilon vs R
    ax8.plot(cond_data['R_um'], cond_data['epsilon'], label=cond, linewidth=1)
    
    # ax9: Temperature (constant for isothermal)
    temp = float(cond.split('=')[1].replace('K', ''))
    ax9.axhline(y=temp, label=cond, linewidth=1)

# Format ax1: R vs time
ax1.set_xlabel(r'Time ($\mu$s)')
ax1.set_ylabel(r'$r_b$ ($\mu$m)')
ax1.set_title('Bubble Radius Evolution')
ax1.set_xscale('log')
ax1.legend(fontsize=7, ncol=2)
ax1.grid(True, alpha=0.3)

# Format ax2: Rdot vs time
ax2.set_xlabel(r'Time ($\mu$s)')
ax2.set_ylabel(r'$\dot{r}_b$ (m/s)')
ax2.set_title('Bubble Wall Velocity')
ax2.set_xscale('log')
ax2.legend(fontsize=7, ncol=2)
ax2.grid(True, alpha=0.3)

# Format ax3: Epsilon vs time
ax3.set_xlabel(r'Time ($\mu$s)')
ax3.set_ylabel(r'$\varepsilon$ (Void Fraction)')
ax3.set_title('Void Fraction Evolution')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.legend(fontsize=7, ncol=2)
ax3.grid(True, alpha=0.3)

# Format ax4: R vs Rdot
ax4.set_xlabel(r'$r_b$ ($\mu$m)')
ax4.set_ylabel(r'$\dot{r}_b$ (m/s)')
ax4.set_title('Bubble Velocity vs Radius')
ax4.legend(fontsize=7, ncol=2)
ax4.grid(True, alpha=0.3)

# Format ax5: rho_m vs time
ax5.set_xlabel(r'Time ($\mu$s)')
ax5.set_ylabel(r'$\rho_m$ (kg/m³)')
ax5.set_title('Mixture Density Evolution')
ax5.set_xscale('log')
ax5.legend(fontsize=7, ncol=2)
ax5.grid(True, alpha=0.3)

# Format ax6: Film thickness vs time
ax6.set_xlabel(r'Time ($\mu$s)')
ax6.set_ylabel(r'Film Thickness ($\mu$m)')
ax6.set_title('Film Thickness Evolution')
ax6.set_xscale('log')
ax6.set_yscale('log')
ax6.legend(fontsize=7, ncol=2)
ax6.grid(True, alpha=0.3)

# Format ax7: P_sat vs time
ax7.set_xlabel(r'Time ($\mu$s)')
ax7.set_ylabel(r'$P_{sat}$ (kPa)')
ax7.set_title('Saturation Pressure (Isothermal)')
ax7.set_xscale('log')
ax7.legend(fontsize=7, ncol=2)
ax7.grid(True, alpha=0.3)

# Format ax8: Epsilon vs R
ax8.set_xlabel(r'$r_b$ ($\mu$m)')
ax8.set_ylabel(r'$\varepsilon$ (Void Fraction)')
ax8.set_title('Void Fraction vs Radius')
ax8.set_yscale('log')
ax8.legend(fontsize=7, ncol=2)
ax8.grid(True, alpha=0.3)

# Format ax9: Temperature
ax9.set_xlabel(r'Time ($\mu$s)')
ax9.set_ylabel(r'Temperature (K)')
ax9.set_title('Droplet Temperature (Isothermal)')
ax9.set_xscale('log')
ax9.legend(fontsize=7, ncol=2)
ax9.grid(True, alpha=0.3)

plt.savefig('song_temp_sweep_c.pdf', dpi=150, bbox_inches='tight')
print("Saved: song_temp_sweep_c.pdf")

# Create additional figure: Rdot vs Film Thickness
fig2, ax = plt.subplots(figsize=(8, 6))
for cond in conditions:
    cond_data = df[df['condition'] == cond]
    ax.plot(cond_data['film_um'], cond_data['Rdot_ms'], label=cond, linewidth=1)
ax.set_xlabel(r'Film Thickness ($\mu$m)')
ax.set_ylabel(r'$\dot{r}_b$ (m/s)')
ax.set_title('Bubble Velocity vs Film Thickness')
ax.set_xlim(50, 0)  # Reversed x-axis
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('song_temp_sweep_c_rdot_vs_film.pdf', dpi=150, bbox_inches='tight')
print("Saved: song_temp_sweep_c_rdot_vs_film.pdf")

# Create additional figure: Rdot vs Time
fig3, ax = plt.subplots(figsize=(8, 6))
for cond in conditions:
    cond_data = df[df['condition'] == cond]
    ax.plot(cond_data['time_us'], cond_data['Rdot_ms'], label=cond, linewidth=1)
ax.set_xlabel(r'Time ($\mu$s)')
ax.set_ylabel(r'$\dot{r}_b$ (m/s)')
ax.set_title('Bubble Wall Velocity Evolution')
ax.set_xscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('song_temp_sweep_c_rdot_vs_time.pdf', dpi=150, bbox_inches='tight')
print("Saved: song_temp_sweep_c_rdot_vs_time.pdf")

# Create additional figure: Radius vs Time
fig4, ax = plt.subplots(figsize=(8, 6))
for cond in conditions:
    cond_data = df[df['condition'] == cond]
    ax.plot(cond_data['time_us'], cond_data['R_um'], label=cond, linewidth=1)
ax.set_xlabel(r'Time ($\mu$s)')
ax.set_ylabel(r'$r_b$ ($\mu$m)')
ax.set_title('Bubble Radius Evolution')
ax.set_xscale('log')
ax.set_xlim(1e-3, 1e2)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.savefig('song_temp_sweep_c_radius_vs_time.pdf', dpi=150, bbox_inches='tight')
print("Saved: song_temp_sweep_c_radius_vs_time.pdf")

# Print summary statistics
print("\n=== Summary Statistics ===")
for cond in conditions:
    cond_data = df[df['condition'] == cond]
    final_time = cond_data['time_us'].iloc[-1]
    final_R = cond_data['R_um'].iloc[-1]
    final_epsilon = cond_data['epsilon'].iloc[-1]
    max_Rdot = cond_data['Rdot_ms'].max()
    print(f"{cond}: t_final={final_time:.2f} μs, R_final={final_R:.2f} μm, "
          f"ε_final={final_epsilon:.4f}, Rdot_max={max_Rdot:.2f} m/s")

print("\nDone!")
