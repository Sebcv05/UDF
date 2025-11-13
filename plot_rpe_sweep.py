#!/usr/bin/env python3
"""
Plot RPE Parameter Sweep Results

Generates two plots:
1. Rdot vs Film Thickness for all T0 values
2. Rdot vs log10(time) for all T0 values

Usage: python3 plot_rpe_sweep.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Use SciencePlots if available (with no-latex to avoid LaTeX errors)
try:
    plt.style.use(['science', 'no-latex', 'ieee'])
    print("Using scienceplots style (no-latex)")
except:
    print("Warning: scienceplots not available, using default style")
    pass

# Read data
print("Loading data...")
df = pd.read_csv('rpe_sweep_all.csv')
df_summary = pd.read_csv('rpe_sweep_summary.csv')

print(f"Loaded {len(df)} data points across {df['T0_K'].nunique()} temperatures")
print(f"\nTemperatures: {sorted(df['T0_K'].unique())} K")

# Get unique temperatures
temps = sorted(df['T0_K'].unique())
n_temps = len(temps)

# Create color map
colors = plt.cm.viridis(np.linspace(0, 1, n_temps))

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Rdot vs Film Thickness
print("\nPlotting Rdot vs Film Thickness...")
for i, T0 in enumerate(temps):
    df_temp = df[df['T0_K'] == T0]
    
    # Get summary data
    summary = df_summary[df_summary['T0_K'] == T0]
    if len(summary) == 0 or summary['max_Rdot_m_s'].values[0] == 0.0:
        print(f"  T0 = {T0:.0f} K: Subcooled - skipping")
        continue
    
    superheat = summary['superheat_K'].values[0]
    
    # Plot the curve (no markers)
    ax1.plot(df_temp['film_thick_um'], abs(df_temp['Rdot_m_s']), 
             color=colors[i], linewidth=2, alpha=0.8,
             label=f'T₀={T0:.0f} K (ΔP={superheat:.2f} bar)')
    
    print(f"  T0 = {T0:.0f} K: plotted")

ax1.set_xlabel('Film Thickness (µm)', fontsize=12, fontweight='bold')
ax1.set_ylabel('|Rdot| (m/s)', fontsize=12, fontweight='bold')
ax1.set_title('Bubble Growth Rate vs Film Thickness\n(P_amb = 2 bar, R₀ = 82.5 µm)', 
              fontsize=13, fontweight='bold')
ax1.legend(fontsize=9, loc='best', framealpha=0.9)
ax1.set_xlim(82.5, 0)  # Reverse x-axis: thick film (left) to thin film (right)
ax1.set_ylim(bottom=0)

# Plot 2: Rdot vs log10(time)
print("\nPlotting Rdot vs log(time)...")
for i, T0 in enumerate(temps):
    df_temp = df[df['T0_K'] == T0]
    
    # Skip subcooled cases
    summary = df_summary[df_summary['T0_K'] == T0]
    if len(summary) == 0 or summary['max_Rdot_m_s'].values[0] == 0.0:
        continue
    
    superheat = summary['superheat_K'].values[0]
    
    # Filter out very early times (numerical artifacts)
    df_temp_filtered = df_temp[df_temp['time_us'] > 1e-6]
    
    # Plot the curve (no markers) - use actual time, not log10
    ax2.plot(df_temp_filtered['time_us'], abs(df_temp_filtered['Rdot_m_s']), 
             color=colors[i], linewidth=2, alpha=0.8,
             label=f'T₀={T0:.0f} K (ΔP={superheat:.2f} bar)')
    
    print(f"  T0 = {T0:.0f} K: plotted")

ax2.set_xlabel('Time (µs)', fontsize=12, fontweight='bold')
ax2.set_ylabel('|Rdot| (m/s)', fontsize=12, fontweight='bold')
ax2.set_title('Bubble Growth Rate vs Time\n(P_amb = 2 bar, R₀ = 82.5 µm)', 
              fontsize=13, fontweight='bold')
ax2.legend(fontsize=9, loc='best', framealpha=0.9)
ax2.set_ylim(bottom=0)
ax2.set_xscale('log')  # Use matplotlib's log scale instead of pre-computed log10

# Add time scale annotations on top axis
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
ax2_top.set_xscale('log')
ax2_top.set_xlabel('Time (µs)', fontsize=10)

plt.tight_layout()

# Save figure
output_file = 'rpe_sweep_plots.pdf'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n=================================================================")
print(f"Plots saved to: {output_file}")
print(f"=================================================================")

plt.show()

# Print summary table
print("\n=================================================================")
print("SUMMARY TABLE")
print("=================================================================")
print(f"{'T0 (K)':<10} {'Superheat':<12} {'Max Rdot':<12} {'Time at Max':<15} {'Film at Max':<15}")
print(f"{'':10} {'(bar)':<12} {'(m/s)':<12} {'(µs)':<15} {'(µm)':<15}")
print("-" * 75)
for _, row in df_summary.iterrows():
    if row['max_Rdot_m_s'] > 0:
        print(f"{row['T0_K']:<10.0f} {row['superheat_K']:<12.3f} "
              f"{row['max_Rdot_m_s']:<12.3f} {row['time_at_max_us']:<15.3e} "
              f"{row['film_thick_at_max_um']:<15.3f}")
    else:
        print(f"{row['T0_K']:<10.0f} {row['superheat_K']:<12.3f} {'SUBCOOLED':<40}")
print("=================================================================")
